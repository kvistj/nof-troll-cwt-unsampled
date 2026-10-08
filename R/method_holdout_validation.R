# Hold-out validation: total CWT recoveries AND PST indicator stock recoveries
#
# Both quantities are computed within the SAME 1000 replicate splits, so the
# total-count and PST-apportioned results are directly comparable and the
# paired advantage comparison is valid across both.
#
# Design
# ------
# Each replicate blinds a random half of the sampled winter strata (and a
# matched random set of summer strata as a control), refits/reassembles both
# estimators on the remaining data, and compares against the known CWT
# recoveries -- both total and PST-only -- actually observed in the blinded
# strata. See prior versions of this script for the full design rationale.
#
# Requires: cwt_fit (with log_exposure), res2, assign_region() defined in
# this session.

library(tidyverse)
library(mgcv)

set.seed(42)

N_REP         <- 1000
SEARCH_WINDOW <- 90
HOLD_FRAC     <- 0.5
KAPPA         <- 50
WINDOW_YRS    <- 0:4

# =====================================================================
# 1. Prepare strata and recovery-level data
# =====================================================================

cwt_fit <- cwt_fit |>
  mutate(season     = ifelse(month %in% c(11,12,1,2,3,4), 'winter', 'summer'),
         date       = as.Date(sprintf('%d-%02d-15', run_year, month)),
         stratum_id = row_number())

rec <- res2 |>
  mutate(age = run_year - brood_year,
         age = pmin(pmax(age, 0), 5)) |>
  filter(age > 1, !is.na(stock)) |>
  mutate(region = assign_region(stock),
         is_pst = region != 'Non-PST stocks') |>
  inner_join(cwt_fit |> select(stratum_id, run_year, month,
                               recovery_location_name, fishery),
             by = c('run_year','month','recovery_location_name','fishery'))

message("Recovery records matched to strata: ", nrow(rec),
        " (", round(100*mean(rec$is_pst), 1), "% PST)")

truth <- rec |>
  group_by(stratum_id) |>
  summarise(n_dec = n(), n_pst = sum(is_pst), .groups = 'drop')

cwt_fit <- cwt_fit |>
  select(-any_of(c('obs_total', 'obs_pst', 'n_dec', 'n_pst'))) |>
  left_join(truth, by = 'stratum_id') |>
  mutate(n_dec     = replace_na(n_dec, 0L),
         n_pst     = replace_na(n_pst, 0L),
         obs_total = cwt_recovs,                      # all tags present
         obs_pst   = case_when(n_dec > 0       ~ cwt_recovs * n_pst / n_dec,
                               cwt_recovs == 0 ~ 0,
                               TRUE            ~ NA_real_))  # tags but none decoded

winter_all <- which(cwt_fit$season == 'winter')
summer_all <- which(cwt_fit$season == 'summer')
n_hold     <- floor(length(winter_all) * HOLD_FRAC)

categories <- rec |> distinct(stock, age) |> arrange(stock, age) |>
  mutate(cat_id = row_number(), region = assign_region(stock),
         is_pst = region != 'Non-PST stocks')
K       <- nrow(categories)
pst_idx <- which(categories$is_pst)

rec <- rec |> left_join(categories |> select(stock, age, cat_id),
                        by = c('stock','age'))

rdirichlet_1 <- function(alpha) {
  g <- rgamma(length(alpha), shape = alpha, rate = 1)
  if (sum(g) == 0) return(alpha / sum(alpha))
  g / sum(g)
}

# =====================================================================
# 2. Tag density model: fit once per replicate, return count AND PST draw
# =====================================================================

predict_td <- function(hold_idx, avail_rec) {
  
  train <- cwt_fit[-hold_idx, ]
  test  <- cwt_fit[ hold_idx, ]
  
  k_month <- max(3, min(8, length(unique(train$month))    - 1))
  k_year  <- max(3, min(8, length(unique(train$run_year)) - 1))
  
  m_h <- tryCatch(
    gam(cwt_recovs ~ s(run_year, k = k_year) + s(month, k = k_month) +
          offset(log_exposure) + recovery_location_name,
        family = tw(link = 'log'), data = train, method = 'REML'),
    error = function(e) NULL)
  if (is.null(m_h)) return(NULL)
  
  cnt <- tryCatch(as.numeric(predict(m_h, newdata = test, type = 'response')),
                  error = function(e) rep(NA_real_, nrow(test)))
  
  pooled <- avail_rec |>
    count(month, recovery_location_name, cat_id, name = 'n') |>
    group_by(month, recovery_location_name) |>
    mutate(p = n / sum(n)) |> ungroup()
  
  pst_pred <- numeric(nrow(test))
  
  for (i in seq_len(nrow(test))) {
    s <- test[i, ]
    
    w <- avail_rec |>
      filter(recovery_location_name == s$recovery_location_name,
             month == s$month,
             run_year %in% (s$run_year - WINDOW_YRS)) |>
      count(cat_id, name = 'n')
    
    pp <- pooled |>
      filter(recovery_location_name == s$recovery_location_name,
             month == s$month)
    
    if (nrow(w) == 0 && nrow(pp) == 0) { pst_pred[i] <- NA_real_; next }
    
    n_vec <- numeric(K); if (nrow(w))  n_vec[w$cat_id]  <- w$n
    p_vec <- numeric(K); if (nrow(pp)) p_vec[pp$cat_id] <- pp$p
    if (sum(p_vec) == 0) p_vec <- n_vec / sum(n_vec)
    
    alpha <- pmax(KAPPA * p_vec, 1e-6) + n_vec
    pst_pred[i] <- cnt[i] * sum(rdirichlet_1(alpha)[pst_idx])
  }
  
  data.frame(stratum_id     = test$stratum_id,
             observed_total = test$obs_total, predicted_total = cnt,
             observed_pst   = test$obs_pst,   predicted_pst   = pst_pred,
             method = 'Tag density model')
}

# =====================================================================
# 3. Ratio estimator: same donor match, return count AND PST prediction
# =====================================================================

predict_ratio <- function(hold_idx, avail_rec) {
  
  test   <- cwt_fit[ hold_idx, ]
  donors <- cwt_fit[-hold_idx, ]
  
  donor_pst <- avail_rec |>
    group_by(stratum_id) |>
    summarise(pst_frac = mean(is_pst), .groups = 'drop')
  donors <- donors |> left_join(donor_pst, by = 'stratum_id')
  
  out_total <- rep(NA_real_, nrow(test))
  out_pst   <- rep(NA_real_, nrow(test))
  
  for (i in seq_len(nrow(test))) {
    s <- test[i, ]
    
    cand <- donors |>
      filter(recovery_location_name == s$recovery_location_name,
             fishery == s$fishery) |>
      mutate(d = as.numeric(abs(date - s$date))) |>
      filter(d <= SEARCH_WINDOW) |>
      arrange(d, date)                              # ties: earlier donor
    
    if (nrow(cand) == 0) next
    
    d1      <- cand[1, ]
    theta_d <- d1$cwt_recovs / d1$number_sampled    # donor tag density
    
    out_total[i] <- s$number_sampled * theta_d      # Eq. 2 at the blinded sample size
    out_pst[i]   <- if (d1$cwt_recovs == 0) 0
    else if (is.na(d1$pst_frac)) NA_real_
    else out_total[i] * d1$pst_frac
  }
  
  data.frame(stratum_id     = test$stratum_id,
             observed_total = test$obs_total, predicted_total = out_total,
             observed_pst   = test$obs_pst,   predicted_pst   = out_pst,
             method = 'Ratio estimator')
}

# =====================================================================
# 4. Replicate splits -- one fit/match per split, both quantities returned
# =====================================================================
run_replicate <- function(rep_id) {
  
  wi <- sample(winter_all, n_hold)
  si <- sample(summer_all, n_hold)
  
  per_set <- function(idx, label) {
    blinded   <- cwt_fit$stratum_id[idx]
    avail_rec <- rec |> filter(!stratum_id %in% blinded)
    
    td <- predict_td(idx, avail_rec)
    re <- predict_ratio(idx, avail_rec)
    if (is.null(td)) return(NULL)
    
    # Score both methods on the same strata: those where the ratio
    # estimator found a donor (mirrors the operational comparison)
    keep <- intersect(re$stratum_id[!is.na(re$predicted_total)],
                      td$stratum_id[!is.na(td$predicted_total)])
    
    res <- bind_rows(td, re) |> filter(stratum_id %in% keep)
    
    keep_pst <- res |>
      group_by(stratum_id) |>
      summarise(ok = all(!is.na(predicted_pst) & !is.na(observed_pst)),
                .groups = 'drop') |>
      filter(ok) |> pull(stratum_id)
    
    res |>
      group_by(method) |>
      summarise(
        n_strata    = n(),
        obs_total   = sum(observed_total),
        pred_total  = sum(predicted_total),
        ratio_total = pred_total / pmax(obs_total, 1),
        
        n_pst     = sum(stratum_id %in% keep_pst),
        obs_pst   = sum(observed_pst[stratum_id %in% keep_pst]),
        pred_pst  = sum(predicted_pst[stratum_id %in% keep_pst]),
        ratio_pst = pred_pst / pmax(obs_pst, 1),
        
        .groups = 'drop'
      ) |>
      mutate(holdout = label, rep = rep_id)
  }
  
  bind_rows(per_set(wi, 'winter'), per_set(si, 'summer'))
}

message("Running ", N_REP, " replicate splits (total + PST, same draws)...")
perf <- bind_rows(lapply(seq_len(N_REP), function(i) {
  if (i %% 50 == 0) message("  replicate ", i, " of ", N_REP)
  run_replicate(i)
}))

# =====================================================================
# 5. Summarise: total and PST, side by side
# =====================================================================

q   <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
rcv <- function(x) IQR(x, na.rm = TRUE) / median(x, na.rm = TRUE)

summarise_quantity <- function(ratio_col, label) {
  perf |>
    group_by(holdout, method) |>
    summarise(
      quantity  = label,
      reps      = n(),
      med_ratio = round(median(.data[[ratio_col]], na.rm = TRUE), 3),
      lwr       = round(q(.data[[ratio_col]], 0.025), 3),
      upr       = round(q(.data[[ratio_col]], 0.975), 3),
      robust_cv = round(rcv(.data[[ratio_col]]), 3),
      .groups = 'drop'
    )
}

summary_tbl <- bind_rows(
  summarise_quantity('ratio_total', 'Total recoveries'),
  summarise_quantity('ratio_pst',   'PST recoveries')
) |> arrange(quantity, holdout, method)

cat("\n=== Calibration across", N_REP, "replicate splits ===\n")
cat("(ratio = summed predicted / summed observed; 1.0 is perfect;\n",
    " robust_cv = IQR / median of the ratio across replicates)\n\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# Paired advantage, both quantities
gap_tbl <- bind_rows(
  perf |> select(rep, holdout, method, ratio = ratio_total) |>
    mutate(quantity = 'Total recoveries'),
  perf |> select(rep, holdout, method, ratio = ratio_pst) |>
    mutate(quantity = 'PST recoveries')
) |>
  pivot_wider(names_from = method, values_from = ratio) |>
  mutate(gap = abs(`Ratio estimator` - 1) - abs(`Tag density model` - 1))

cat("\n=== Advantage of tag density model (positive favours it) ===\n")
print(as.data.frame(
  gap_tbl |> group_by(quantity, holdout) |>
    summarise(median_advantage = round(median(gap, na.rm = TRUE), 3),
              pct_favouring    = round(100*mean(gap > 0, na.rm = TRUE), 1),
              .groups = 'drop')), row.names = FALSE)

# =====================================================================
# 6. Plot: 4 panels -- quantity (total/PST) x holdout (summer/winter)
# =====================================================================
plot_dat <- bind_rows(
  perf |> transmute(rep, holdout, method, ratio = ratio_total,
                    quantity = 'Total recoveries'),
  perf |> transmute(rep, holdout, method, ratio = ratio_pst / ratio_total,
                    quantity = 'PST share'),
  perf |> transmute(rep, holdout, method, ratio = ratio_pst,
                    quantity = 'PST recoveries')
) |>
  mutate(holdout  = recode(holdout, winter = 'Winter hold-out (Nov–Apr)',
                           summer = 'Summer hold-out (control)'),
         quantity = factor(quantity, levels = c('Total recoveries',
                                                'PST share',
                                                'PST recoveries')),
         holdout  = factor(holdout, levels = c('Summer hold-out (control)',
                                               'Winter hold-out (Nov–Apr)')))

p <- ggplot(plot_dat, aes(x = ratio, fill = method)) +
  geom_histogram(bins = 40, alpha = 0.6, position = 'identity',
                 color = 'black', linewidth = 0.2) +
  geom_vline(xintercept = 1, lty = 2, color = 'grey20') +   # drawn on top
  facet_grid(quantity ~ holdout, scales = 'fixed') +
  scale_fill_manual(values = c('Ratio estimator'   = '#440154',
                               'Tag density model' = '#E5B84B')) +
  labs(x = 'Predicted / observed', y = 'Replicate splits', fill = 'Method') +
  theme_bw() +
  coord_cartesian(xlim = c(0, 3)) +
  theme(text = element_text(size = 14), legend.position = 'bottom')

print(p)
ggsave("out/holdout_validation_combined.png", p, height = 7, width = 7, dpi = 300)


perf |>
  mutate(share   = ratio_pst / ratio_total,
         err_log = abs(log(share)),
         err_abs = abs(share - 1)) |>
  select(rep, holdout, method, err_log, err_abs) |>
  pivot_wider(names_from = method, values_from = c(err_log, err_abs)) |>
  group_by(holdout) |>
  summarise(
    pct_closer_log = 100 * mean(`err_log_Tag density model` < `err_log_Ratio estimator`, na.rm = TRUE),
    pct_closer_abs = 100 * mean(`err_abs_Tag density model` < `err_abs_Ratio estimator`, na.rm = TRUE),
    med_share_tdm  = median(perf$ratio_pst[perf$method == "Tag density model" & perf$holdout == first(holdout)] /
                              perf$ratio_total[perf$method == "Tag density model" & perf$holdout == first(holdout)], na.rm = TRUE),
    med_share_re   = median(perf$ratio_pst[perf$method == "Ratio estimator" & perf$holdout == first(holdout)] /
                              perf$ratio_total[perf$method == "Ratio estimator" & perf$holdout == first(holdout)], na.rm = TRUE)
  )
