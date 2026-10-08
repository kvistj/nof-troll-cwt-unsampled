# Prior sensitivity for the Dirichlet composition layer (Table S1, Figure 9)
#
# Replaces method2_prior_comparison.R. Table S1 previously came from an older
# run that refit the GAM by bootstrap, so its kappa = 50 row did not match the
# headline estimate from posterior simulation. This script uses the same
# posterior rate-model draws as method2_uncertainty.R and reproduces the
# headline kappa = 50 result exactly.
#
# How the match is guaranteed:
#   1. The rate-model draws are generated with the same seed and the same
#      calls as draw_counts_posterior(), so count_draws is identical.
#   2. The RNG state immediately after those draws is saved and restored
#      before each prior, so every prior sees the same Dirichlet random
#      stream the headline run saw. Differences among rows are therefore
#      due to the prior alone, not Monte Carlo noise.
#
# Usage (fresh session, after running the .qmd through the composition chunks
# so that cwt_fit, m, pred_dat and res2 exist):
#   source('method2_prior_sensitivity.R')

library(tidyverse)
library(mgcv)

# ---------------------------------------------------------------------
# 0. Build the shared objects via method2_uncertainty.R (posterior only)
# ---------------------------------------------------------------------
# Settings are set explicitly here so stale values from an earlier session
# (e.g. kappa left at 200 from a sweep) cannot leak in.
n_draws    <- 1000
kappa      <- 50
gam_method <- 'posterior'
source('R/method2_uncertainty.R')

needed <- c('pred_strata', 'n_mat', 'pooled_mat', 'has_comp', 'pst_idx',
            'K', 'n_strata', 'rdirichlet_1', 'm', 'results')
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing)) stop('Missing objects: ', paste(missing, collapse = ', '))

headline <- results$posterior          # the run behind the reported 1,158

# ---------------------------------------------------------------------
# 1. Rate-model draws: identical to draw_counts_posterior()
# ---------------------------------------------------------------------
set.seed(42)
count_draws      <- count_draws_post
rng_after_counts <- rng_before_dirichlet      # restored before every prior

# ---------------------------------------------------------------------
# 2. Priors
# ---------------------------------------------------------------------
# Posterior alpha is built exactly as in method2_uncertainty.R.
make_alpha <- function(prior_mat) {
  a <- n_mat + prior_mat
  a[a <= 0] <- 1e-6
  a
}

kappas <- c(5, 10, 25, 50, 100, 200)

priors <- c(
  setNames(lapply(kappas, function(k) make_alpha(k * pooled_mat)),
           paste0('eb_', kappas)),
  list(jeffreys = make_alpha(matrix(0.5, n_strata, K)),
       uniform  = make_alpha(matrix(1.0, n_strata, K)))
)

prior_info <- tibble(
  key   = names(priors),
  prior = c(rep('Empirical Bayes', length(kappas)),
            'Jeffreys (αk = 1/2)', 'Uniform (αk = 1)'),
  kappa = c(kappas, NA, NA)
)

# ---------------------------------------------------------------------
# 3. Apportion under each prior with a common random stream
# ---------------------------------------------------------------------
apportion_pst <- function(alpha_post) {
  pst <- numeric(n_draws)
  for (b in seq_len(n_draws)) {
    cw <- count_draws[b, ]
    s  <- 0
    for (i in seq_len(n_strata)) {
      if (!has_comp[i] || is.na(cw[i])) next
      pd <- rdirichlet_1(alpha_post[i, ])
      s  <- s + cw[i] * sum(pd[pst_idx])
    }
    pst[b] <- s
  }
  pst
}

pst_by_prior <- lapply(names(priors), function(k) {
  message('  ', k)
  assign('.Random.seed', rng_after_counts, envir = .GlobalEnv)
  apportion_pst(priors[[k]])
})
names(pst_by_prior) <- names(priors)

# ---------------------------------------------------------------------
# 4. Check: kappa = 50 must reproduce the headline run exactly
# ---------------------------------------------------------------------
stopifnot(isTRUE(all.equal(pst_by_prior$eb_50, headline$total_pst)))
message('kappa = 50 row reproduces the headline PST distribution exactly.')

# ---------------------------------------------------------------------
# 5. Table S1
# ---------------------------------------------------------------------
fmt <- function(x) formatC(round(x), format = 'd', big.mark = ',')

table_s1 <- prior_info |>
  mutate(
    lower  = map_dbl(key, ~ quantile(pst_by_prior[[.x]], 0.025)),
    median = map_dbl(key, ~ quantile(pst_by_prior[[.x]], 0.500)),
    upper  = map_dbl(key, ~ quantile(pst_by_prior[[.x]], 0.975)),
    pct_vs_k50 = 100 * (median / median[key == 'eb_50'] - 1)
  )

print(table_s1 |> select(prior, kappa, median, lower, upper, pct_vs_k50),
      n = Inf, digits = 5)

table_s1_out <- table_s1 |>
  transmute(Prior = prior,
            `κ` = ifelse(is.na(kappa), '–', as.character(kappa)),
            Median = fmt(median),
            `95% CI` = paste0(fmt(lower), '–', fmt(upper)))

# write.csv(table_s1_out, 'out/table_s1_prior_sensitivity.csv', row.names = FALSE)

# Values quoted in the text
eb_range <- range(table_s1$pct_vs_k50[!is.na(table_s1$kappa)])
flat_pct <- table_s1$pct_vs_k50[is.na(table_s1$kappa)]
cat(sprintf('\nEmpirical Bayes medians span %.1f%% to %.1f%% of kappa = 50\n',
            eb_range[1], eb_range[2]))
cat(sprintf('Jeffreys / uniform medians are %.0f%% / %.0f%% higher than kappa = 50\n',
            flat_pct[1], flat_pct[2]))

# ---------------------------------------------------------------------
# 6. Figure 9
# ---------------------------------------------------------------------
total_all <- headline$total_cwt
q_all     <- quantile(total_all, c(0.025, 0.5, 0.975))

p_all <- ggplot(data.frame(total = total_all), aes(x = total)) +
  geom_histogram(fill = 'azure2', color = 'black', bins = 40) +
  geom_vline(xintercept = q_all[2], lty = 2, linewidth = 0.8, alpha = 0.6) +
  geom_vline(xintercept = q_all[c(1, 3)], lty = 3, linewidth = 0.6, alpha = 0.6) +
  coord_cartesian(xlim = c(700, 1800)) +
  labs(x = '', y = 'Count',
       title = 'All stocks') +
  theme_bw() + theme(text = element_text(size = 14))

fig_keys <- c('eb_50', 'jeffreys', 'uniform')
fig_labs <- c(eb_50 = 'Empirical Bayes (κ = 50) ',
              jeffreys = 'Jeffreys (αk = 1/2)',
              uniform  = 'Uniform (αk = 1)')

plot_dat <- bind_rows(lapply(fig_keys, function(k)
  data.frame(prior = fig_labs[[k]], total = pst_by_prior[[k]]))) |>
  mutate(prior = factor(prior, levels = fig_labs[fig_keys]))

med_dat <- plot_dat |> group_by(prior) |> summarise(med = median(total))

p_pst <- ggplot(plot_dat, aes(x = total, fill = prior)) +
  geom_histogram(color = 'black', bins = 40, alpha = 0.6,
                 position = 'identity') +
  geom_vline(data = med_dat, aes(xintercept = med, color = prior),
             lty = 2, linewidth = 0.8, show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set2') +
  scale_color_brewer(palette = 'Set2') +
  coord_cartesian(xlim = c(700, 1800)) +
  labs(x = 'Estimated CWT recoveries', y = 'Count',
       fill = 'Dirichlet prior', title = 'PST indicator stocks') +
  theme_bw() + theme(text = element_text(size = 14), legend.position = 'bottom')

cowplot::plot_grid(p_all, p_pst, ncol = 1, labels = 'auto',
                   rel_heights = c(1, 1.15))

ggsave('out/fig9_prior_comparison.png', height = 6, width = 8, dpi = 300)

