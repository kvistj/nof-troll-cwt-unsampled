# Method 2 Uncertainty Estimation
#
# Propagates two sources of uncertainty to the estimate of total CWT recoveries
# from unsampled strata:
#
#   (1) Rate model  -- uncertainty in the fitted tag density
#   (2) Composition -- Dirichlet posterior over stock/age composition
#
# The rate-model contribution can be obtained two ways, controlled by
# `gam_method` below:
#
#   'posterior'  Draw coefficients from the approximate posterior implied by
#                the fitted GAM, beta | data ~ N(beta_hat, V_beta), and push
#                them through the linear predictor. No refitting. Combined with
#                the Dirichlet layer this yields credible intervals throughout.
#
#   'bootstrap'  Resample the training strata with replacement and refit the
#                GAM each iteration. Captures influential-stratum effects and
#                some model misspecification that posterior simulation does not.
#
# Running both and comparing is informative: close agreement supports reporting
# the posterior-simulation intervals as credible intervals; material divergence
# indicates the resampling variability is not purely parametric and warrants
# investigation before choosing.
#
# Usage: source after running the .qmd data prep through the composition
# chunks, so that cwt_fit, m, pred_dat, and res2 are in the environment.

library(tidyverse)
library(mgcv)

# --- Configuration ---
n_draws    <- 1000
kappa      <- 50          # empirical-Bayes prior concentration
if (!exists('gam_method'))
  gam_method <- 'both'     # 'posterior', 'bootstrap', or 'both'

set.seed(42)

# =====================================================================
# 1. Fixed stock/age category grid
# =====================================================================
# Every stratum's count vector spans the SAME category set, with structural
# zeros where a category was not observed, so that the prior can place mass on
# categories absent from a given stratum.

res_aged <- res2 |>
  mutate(age = run_year - brood_year,
         age = ifelse(age < 0, 0, age),
         age = ifelse(age >= 5, 5, age)) |>
  filter(age > 1, !is.na(stock))

stopifnot(all(c('run_year','month','recovery_location_name',
                'stock','brood_year') %in% names(res2)))

categories <- res_aged |>
  distinct(stock, age) |>
  arrange(stock, age) |>
  mutate(cat_id = row_number())

K <- nrow(categories)
message("Category grid: K = ", K, " stock/age combinations")

assign_region <- function(stock) {
  case_when(
    stock %in% c('ELW','GAD','NIS','NSF','GRO','GRN','HOK','SAM','SKF','SSF',
                 'STL','SKY','WHY','WRY','WHF','WRF','SKS','SPY') ~ "Puget Sound",
    stock %in% c('QNT','QN2','SOO','QUE','QUI')                   ~ 'WA Coast',
    stock %in% c('SUM','URB','SPR','LYY','LYF','HAN','KLM','MSH',
                 'NIC','SHU','SMK')                               ~ 'Columbia River',
    stock %in% c('ATN','BQR','CHI','CHU','COW','CKO','HAR','RBT',
                 'NAN','DOM')                                     ~ 'Canada',
    stock %in% c('LRH','WSH','CWF','LRW','ELK','SRH','SRN')       ~ 'Oregon',
    stock %in% c('ALP')                                           ~ 'Alaska',
    TRUE                                                          ~ 'Non-PST stocks'
  )
}

categories$region <- assign_region(categories$stock)
pst_idx <- which(categories$region != 'Non-PST stocks')

# =====================================================================
# 2. Windowed and pooled COUNT tables
# =====================================================================
# Lagging five-year window: recoveries from years t-4 through t inform
# composition for window-year t.

window_counts <- res_aged |>
  cross_join(tibble(year_offset = 0:4)) |>
  mutate(window_year = run_year + year_offset) |>
  count(window_year, month, recovery_location_name, stock, age, name = 'n')

pooled_counts <- res_aged |>
  count(month, recovery_location_name, stock, age, name = 'n')

pooled_props <- pooled_counts |>
  group_by(month, recovery_location_name) |>
  mutate(pooled_p = n / sum(n)) |>
  ungroup() |>
  select(month, recovery_location_name, stock, age, pooled_p)

# =====================================================================
# 3. Per-stratum count and prior matrices
# =====================================================================

pred_strata <- pred_dat |>
  filter(type == 'Predicted') |>
  mutate(stratum_id = row_number()) |>
  select(stratum_id, run_year, month, recovery_location_name,
         fishery, number_caught, log_exposure) |>
  as.data.frame()

n_strata <- nrow(pred_strata)

n_mat      <- matrix(0, nrow = n_strata, ncol = K)
pooled_mat <- matrix(0, nrow = n_strata, ncol = K)
has_comp   <- logical(n_strata)

cat_key <- categories |> select(stock, age, cat_id)

for (i in seq_len(n_strata)) {
  s <- pred_strata[i, ]
  
  w <- window_counts |>
    filter(window_year == s$run_year, month == s$month,
           recovery_location_name == s$recovery_location_name)
  
  if (nrow(w) == 0) {
    w <- pooled_counts |>
      filter(month == s$month,
             recovery_location_name == s$recovery_location_name)
  }
  if (nrow(w) == 0) { has_comp[i] <- FALSE; next }
  has_comp[i] <- TRUE
  
  w <- left_join(w, cat_key, by = c('stock','age'))
  n_mat[i, w$cat_id] <- w$n
  
  pp <- pooled_props |>
    filter(month == s$month,
           recovery_location_name == s$recovery_location_name) |>
    left_join(cat_key, by = c('stock','age'))
  
  if (nrow(pp) > 0) {
    pooled_mat[i, pp$cat_id] <- pp$pooled_p
  } else {
    pooled_mat[i, ] <- n_mat[i, ] / sum(n_mat[i, ])
  }
}

message("Predicted strata: ", n_strata,
        " | with composition: ", sum(has_comp),
        " | without: ", sum(!has_comp))

# =====================================================================
# 4. Samplers
# =====================================================================

rdirichlet_1 <- function(alpha) {
  g <- rgamma(length(alpha), shape = alpha, rate = 1)
  if (sum(g) == 0) return(alpha / sum(alpha))
  g / sum(g)
}

# Posterior alpha, fixed across draws
alpha_post <- sweep(n_mat, 2, 0, '+') + kappa * pooled_mat
alpha_post[alpha_post <= 0] <- 1e-6

# ---------------------------------------------------------------------
# Rate-model draws: posterior simulation
# ---------------------------------------------------------------------
# beta | data ~ N(beta_hat, V_beta). Using unconditional = TRUE incorporates
# uncertainty in the smoothing parameters. lpmatrix excludes the offset, which
# is added on the link scale before exponentiating.
draw_counts_posterior <- function(n_draws) {
  Xp <- predict(m, newdata = pred_strata, type = 'lpmatrix')
  Vb <- vcov(m, unconditional = TRUE)
  L  <- chol(Vb)
  bd <- matrix(rnorm(n_draws * ncol(Vb)), n_draws) %*% L
  bd <- sweep(bd, 2, coef(m), '+')
  eta <- bd %*% t(Xp)
  exp(sweep(eta, 2, pred_strata$log_exposure, '+'))
}

# ---------------------------------------------------------------------
# Rate-model draws: nonparametric bootstrap
# ---------------------------------------------------------------------
draw_counts_bootstrap <- function(n_draws) {
  out <- matrix(NA_real_, nrow = n_draws, ncol = n_strata)
  fails <- 0
  for (b in seq_len(n_draws)) {
    if (b %% 100 == 0) message("    refit ", b, " of ", n_draws)
    ok <- tryCatch({
      bd <- cwt_fit[sample(nrow(cwt_fit), replace = TRUE), ]
      mb <- gam(cwt_recovs ~ s(run_year) + s(month) + offset(log_exposure) +
                  recovery_location_name,
                family = tw(link = 'log'), data = bd, method = 'REML')
      out[b, ] <- as.numeric(predict(mb, newdata = pred_strata,
                                     type = 'response'))
      TRUE
    }, error = function(e) FALSE)
    if (!ok) fails <- fails + 1
  }
  attr(out, 'failures') <- fails
  out
}

# =====================================================================
# 5. Apportionment
# =====================================================================
# Given a matrix of rate-model draws (rows = draws, cols = strata), draw a
# composition per stratum per row and return totals.

apportion <- function(count_draws) {
  nd <- nrow(count_draws)
  tot_all <- numeric(nd); tot_pst <- numeric(nd); tot_app <- numeric(nd)
  
  for (b in seq_len(nd)) {
    cw <- count_draws[b, ]
    if (all(is.na(cw))) { tot_all[b] <- NA; tot_pst[b] <- NA; next }
    tot_all[b] <- sum(cw, na.rm = TRUE)
    ps <- 0; ap <- 0
    for (i in seq_len(n_strata)) {
      if (!has_comp[i] || is.na(cw[i])) next
      pd <- rdirichlet_1(alpha_post[i, ])
      ps <- ps + cw[i] * sum(pd[pst_idx])
      ap <- ap + cw[i]
    }
    tot_pst[b] <- ps; tot_app[b] <- ap
  }
  data.frame(total_cwt = tot_all, total_apport = tot_app, total_pst = tot_pst)
}

# =====================================================================
# 6. Run
# =====================================================================

results <- list()

if (gam_method %in% c('posterior','both')) {
  message("Posterior simulation (", n_draws, " draws)...")
  t0 <- Sys.time()
  count_draws_post     <- draw_counts_posterior(n_draws)
  rng_before_dirichlet <- .Random.seed   # reused by the prior sensitivity script
  results$posterior    <- apportion(count_draws_post)
  message("  done in ", round(difftime(Sys.time(), t0, units = 'secs'), 1), "s")
}

if (gam_method %in% c('bootstrap','both')) {
  set.seed(4242)   # bootstrap stream independent of whether the posterior ran
  message("Nonparametric bootstrap (", n_draws, " refits)...")
  t0 <- Sys.time()
  cd <- draw_counts_bootstrap(n_draws)
  message("  ", attr(cd, 'failures'), " iteration(s) failed to converge")
  results$bootstrap <- apportion(cd)
  message("  done in ", round(difftime(Sys.time(), t0, units = 'mins'), 1), " min")
}
# =====================================================================
# 7. Compare
# =====================================================================

ci <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

summarise_run <- function(df, label) {
  a <- ci(df$total_cwt); p <- ci(df$total_pst)
  data.frame(
    method   = label,
    quantity = c('All stocks','PST indicator stocks'),
    lower    = c(a[1], p[1]),
    median   = c(a[2], p[2]),
    upper    = c(a[3], p[3]),
    width    = c(a[3]-a[1], p[3]-p[1])
  )
}

cmp <- bind_rows(lapply(names(results),
                        function(k) summarise_run(results[[k]], k)))

cat("\n=== Point estimate ===\n")
cat("All stocks:", sum(pred_dat$cwt_recovs[pred_dat$type == 'Predicted'],
                       na.rm = TRUE), "\n")

cat("\n=== Interval comparison ===\n")
print(cmp, row.names = FALSE, digits = 5)

if (length(results) == 2) {
  wp <- cmp$width[cmp$method == 'posterior']
  wb <- cmp$width[cmp$method == 'bootstrap']
  cat("\nWidth ratio (bootstrap / posterior):\n")
  cat("  All stocks: ", round(wb[1]/wp[1], 3), "\n")
  cat("  PST stocks: ", round(wb[2]/wp[2], 3), "\n")
  cat("\nRatios near 1 indicate the resampling variability is essentially\n",
      "parametric, supporting the posterior-simulation intervals. Ratios\n",
      "materially above 1 suggest influential strata or misspecification\n",
      "that posterior simulation does not capture.\n", sep = "")
}

# =====================================================================
# 8. Plot
# =====================================================================

plot_dat <- bind_rows(lapply(names(results), function(k)
  results[[k]] |> mutate(method = k))) |>
  pivot_longer(c(total_cwt, total_pst),
               names_to = 'quantity', values_to = 'total') |>
  mutate(quantity = recode(quantity,
                           total_cwt = 'All stocks',
                           total_pst = 'PST indicator stocks'))

# Reference lines from the posterior run: the reported credible intervals
vlines <- plot_dat |>
  filter(method == 'posterior') |>
  group_by(quantity) |>
  summarise(Median   = median(total, na.rm = TRUE),
            lo       = quantile(total, 0.025, na.rm = TRUE),
            hi       = quantile(total, 0.975, na.rm = TRUE),
            .groups  = 'drop') |>
  pivot_longer(c(Median, lo, hi), names_to = 'stat', values_to = 'x') |>
  mutate(stat = ifelse(stat == 'Median', 'Median', '95% CI'))

mk_panel <- function(q, show_x = TRUE) {
  d <- filter(plot_dat, quantity == q)
  v <- filter(vlines,   quantity == q)
  ggplot(d, aes(x = total, fill = method)) +
    geom_histogram(color = 'black', bins = 40, alpha = 0.55,
                   position = 'identity') +
    geom_vline(data = v, aes(xintercept = x, linetype = stat),
               color = 'grey20', linewidth = 0.5) +
    scale_fill_manual(values = c("bootstrap" = "#440154", "posterior" = "#E5B84B")) +
    scale_linetype_manual(values = c('Median' = 'dashed', '95% CI' = 'dotted')) +
    coord_cartesian(xlim = range(plot_dat$total, na.rm = TRUE)) +
    labs(x = if (show_x) 'Total estimated CWT recoveries' else NULL,
         y = 'Count', title = q,
         fill = 'Uncertainty method', linetype = NULL) +
    theme_bw() +
    theme(text = element_text(size = 15),
          legend.position = 'bottom',
          legend.box = 'horizontal')
}

p_all <- mk_panel('All stocks', show_x = FALSE)
p_pst <- mk_panel('PST indicator stocks', show_x = TRUE)

# Pull the legend from one panel, then strip legends from both
legend <- cowplot::get_legend(p_all)

panels <- cowplot::plot_grid(
  p_all + theme(legend.position = 'none'),
  p_pst + theme(legend.position = 'none'),
  ncol = 1, labels = 'auto', align = 'v'
)

cowplot::plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.1))

if (length(results) == 2) {
  cowplot::plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.1))
  ggsave("out/method2_uncertainty_comparison.png", height = 7.5, width = 8, dpi = 300)
}