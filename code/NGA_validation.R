# =============================================================================
# NGA_validation.R
# GATING PHASE — must pass before NGA_run_highrisk.R is trusted.
# "Show the model is good first": validate the joint coregionalized GP hurdle
# model extensively, with special attention to the TAIL (the >20,000 question
# is a tail extrapolation, so tail adequacy is decisive).
#
# Sections:
#   A. Convergence diagnostics (Rhat, ESS, divergences, E-BFMI, treedepth)
#   B. Prevalence: spatial block CV -> AUC, Brier, log-loss, reliability
#   C. Intensity: spatial-block CV PIT, CRPS, interval coverage
#   D. Tail-specific: GPD mean-residual-life, shape stability, PP of the maximum
#   E. Residual spatial autocorrelation: Moran's I + empirical variogram
#   F. Model comparison: LOO / WAIC across tails x {indep vs coregion} x {+/-cov}
#
# Outputs figures + metrics to results/highrisk_20000/validation/.
# =============================================================================

source("code/NGA_common.R")
suppressMessages({ library(loo); library(pROC); library(scoringRules); library(ggplot2) })

VAL_DIR <- file.path(CFG$paths$out_dir, "validation")
dir.create(VAL_DIR, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")
metrics <- list()

dat <- load_nga()
models <- compile_models()

# Primary body model used throughout validation (swap to compare).
PRIMARY <- "lognormal"
log_msg("Fitting primary model for validation:", PRIMARY)
fit <- fit_hurdle(dat, PRIMARY, models$hurdle)
saveRDS(fit, file.path(VAL_DIR, "fit_primary.rds"))

# -----------------------------------------------------------------------------
# A. Convergence
# -----------------------------------------------------------------------------
log_msg("A. Convergence diagnostics...")
sm <- rstan::summary(fit)$summary
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
metrics$convergence <- tibble(
  max_rhat       = max(sm[, "Rhat"], na.rm = TRUE),
  min_bulk_ess   = min(sm[, "n_eff"], na.rm = TRUE),
  n_divergent    = sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
  max_treedepth_hit = sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= CFG$max_treedepth))),
  ebfmi_min      = min(rstan::get_bfmi(fit))
)
print(metrics$convergence)

# -----------------------------------------------------------------------------
# Spatial block CV scaffold: K folds by k-means on village coordinates.
# Refit on training villages, predict held-out village individuals via kriging.
# -----------------------------------------------------------------------------
make_folds <- function(villages, K = 10, seed = 1) {
  set.seed(seed)
  km <- kmeans(as.matrix(villages[, c("utmx", "utmy")]), centers = K, nstart = 5)
  villages$fold <- km$cluster
  villages
}

predict_holdout <- function(fit_tr, dat_tr, hold_coords, X_hold, family) {
  # ICM kriging (conditional mean) of f_p / f_c to held-out villages. The shared
  # correlation K(rho) makes each field's conditional mean a unit-variance krige
  # (the marginal SD cancels in the mean), then add intercept + covariates.
  post  <- rstan::extract(fit_tr)
  Dgv_h <- fields_rdist(hold_coords, as.matrix(dat_tr$villages[, c("utmx", "utmy")]))
  f_p_h <- krige_mean(post$f_p, dat_tr$Dmat, Dgv_h, post$rho)
  f_c_h <- krige_mean(post$f_c, dat_tr$Dmat, Dgv_h, post$rho)
  eta_p_h <- sweep(f_p_h, 1, post$b0_p, "+")
  eta_c_h <- sweep(f_c_h, 1, post$b0_c, "+")
  if (ncol(X_hold) > 0) {
    eta_p_h <- eta_p_h + post$beta_p %*% t(X_hold)   # (S x P)(P x H) = S x H
    eta_c_h <- eta_c_h + post$beta_c %*% t(X_hold)
  }
  list(pi = plogis(eta_p_h), eta_c = eta_c_h,
       shape = post$shape, sigma_ln = post$sigma_ln)
}

run_spatial_cv <- function(dat, family, K = 6, iter = 800, warmup = 400) {
  vill <- make_folds(dat$villages, K)
  preds <- list()
  for (k in seq_len(K)) {
    log_msg("  CV fold", k, "/", K)
    tr_vill <- vill |> filter(fold != k)
    # rebuild a training data object restricted to training villages
    dtr <- dat$d |> filter(vill_key %in% tr_vill$vill_key)
    tr_vill <- tr_vill |> mutate(vid = row_number())
    dtr <- dtr |> select(-vid) |> left_join(tr_vill |> select(vill_key, vid), by = "vill_key")
    Xtr <- if (ncol(dat$Xv) > 0) dat$Xv[dtr$vill_key, , drop = FALSE] else matrix(0, nrow(dtr), 0)
    dat_tr <- list(d = dtr, villages = tr_vill, V = nrow(tr_vill),
                   Dmat = as.matrix(dist(as.matrix(tr_vill[, c("utmx","utmy")]))),
                   X = Xtr, pos_idx = which(dtr$y_binary == 1))
    fit_tr <- fit_hurdle(dat_tr, family, models$hurdle, iter = iter, warmup = warmup)
    ho_vill <- vill |> filter(fold == k)
    ho_ind  <- dat$d |> filter(vill_key %in% ho_vill$vill_key)
    Xho <- if (ncol(dat$Xv) > 0) dat$Xv[ho_vill$vill_key, , drop = FALSE] else matrix(0, nrow(ho_vill), 0)
    p <- predict_holdout(fit_tr, dat_tr,
                         as.matrix(ho_vill[, c("utmx","utmy")]), Xho, family)
    # map held-out individuals to held-out village row
    vmap <- match(ho_ind$vill_key, ho_vill$vill_key)
    preds[[k]] <- list(ind = ho_ind, vmap = vmap, p = p)
  }
  preds
}

# -----------------------------------------------------------------------------
# B. Prevalence CV metrics
# -----------------------------------------------------------------------------
log_msg("B. Spatial CV for prevalence (this refits K models)...")
cv <- run_spatial_cv(dat, PRIMARY, K = 6)

yb <- unlist(lapply(cv, function(f) f$ind$y_binary))
pi_hat <- unlist(lapply(cv, function(f) colMeans(f$p$pi)[f$vmap]))
roc_obj <- pROC::roc(yb, pi_hat, quiet = TRUE)
metrics$prevalence_cv <- tibble(
  auc    = as.numeric(pROC::auc(roc_obj)),
  brier  = mean((pi_hat - yb)^2),
  logloss= -mean(yb * log(pi_hat) + (1 - yb) * log(1 - pi_hat))
)
print(metrics$prevalence_cv)

# Reliability diagram
rel <- tibble(p = pi_hat, y = yb) |>
  mutate(bin = cut(p, seq(0, 1, 0.1), include.lowest = TRUE)) |>
  group_by(bin) |> summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop")
ggplot(rel, aes(pred, obs, size = n)) +
  geom_abline(linetype = 2) + geom_point() +
  labs(title = "Prevalence reliability (spatial CV)", x = "Predicted", y = "Observed") +
  theme_minimal()
ggsave(file.path(VAL_DIR, "B_reliability.png"), width = 5, height = 5, dpi = 150)

# -----------------------------------------------------------------------------
# C. Intensity CV metrics: CRPS, coverage, spatial-block CV PIT
# -----------------------------------------------------------------------------
log_msg("C. Intensity calibration (CRPS, coverage, spatial-block CV PIT)...")
# Predictive samples for held-out positives under the family
intensity_metrics <- map_dfr(cv, function(f) {
  pos <- which(f$ind$y_binary == 1)
  if (!length(pos)) return(tibble())
  map_dfr(pos, function(j) {
    v <- f$vmap[j]
    eta <- f$p$eta_c[, v]
    draws <- if (PRIMARY == "lognormal") rlnorm(length(eta), eta, f$p$sigma_ln)
             else if (PRIMARY == "weibull") rweibull(length(eta), f$p$shape, exp(eta))
             else rgamma(length(eta), shape = f$p$shape, rate = f$p$shape / exp(eta))
    obs <- f$ind$y_intensity[j]
    tibble(obs = obs,
           crps = scoringRules::crps_sample(obs, draws),
           pit  = mean(draws <= obs),
           in90 = obs >= quantile(draws, .05) & obs <= quantile(draws, .95))
  })
})
metrics$intensity_cv <- tibble(
  mean_crps = mean(intensity_metrics$crps),
  cover90   = mean(intensity_metrics$in90)
)
print(metrics$intensity_cv)
write_csv(intensity_metrics, file.path(VAL_DIR, "intensity_cv_metrics.csv"))
pit_plot <- ggplot(intensity_metrics, aes(pit)) +
  geom_histogram(boundary = 0, bins = 10, fill = "grey70", colour = "white") +
  geom_hline(yintercept = nrow(intensity_metrics) / 10, linetype = 2) +
  labs(x = "PIT", y = "Count") +
  theme_minimal(base_size = 11)
# Preserve the legacy filename but export publication-quality raster and vectors.
ggsave(file.path(VAL_DIR, "C_loo_pit.png"), pit_plot, width = 5, height = 4, dpi = 600)
ggsave(file.path(VAL_DIR, "C_loo_pit.pdf"), pit_plot, width = 5, height = 4)
ggsave(file.path(VAL_DIR, "C_loo_pit.eps"), pit_plot, width = 5, height = 4,
       device = grDevices::cairo_ps)

# -----------------------------------------------------------------------------
# D. Tail-specific validation (decisive for the 20,000 question)
# -----------------------------------------------------------------------------
log_msg("D. Tail diagnostics...")
yc <- dat$d$y_intensity[dat$pos_idx]

# D1. Mean residual life plot (should be ~linear above a valid POT threshold u)
mrl <- map_dfr(seq(100, quantile(yc, .98), length.out = 40), function(u) {
  ex <- yc[yc > u] - u
  tibble(u = u, mrl = mean(ex), se = sd(ex) / sqrt(length(ex)), n = length(ex))
})
ggplot(mrl, aes(u, mrl)) + geom_line() +
  geom_ribbon(aes(ymin = mrl - 2*se, ymax = mrl + 2*se), alpha = .2) +
  geom_vline(xintercept = CFG$gpd_u, linetype = 2, colour = "red") +
  labs(title = "Mean residual life (GPD threshold choice)", x = "u (mf/mL)", y = "mean excess") +
  theme_minimal()
ggsave(file.path(VAL_DIR, "D_mrl.png"), width = 6, height = 4, dpi = 150)

# D2. GPD shape stability across thresholds (xi should be ~constant above valid u)
fit_gpd_at <- function(u) {
  e <- yc[yc > u] - u
  if (length(e) < 10) return(c(u = u, xi = NA, sigma = NA))
  f <- fit_gpd(list(d = tibble(y_intensity = yc), pos_idx = seq_along(yc)),
               models$gpd, modifyList(CFG, list(gpd_u = u, iter = 1500, warmup = 750)))
  c(u = u, xi = mean(rstan::extract(f)$xi), sigma = mean(rstan::extract(f)$sigma))
}
shape_stab <- as_tibble(do.call(rbind, lapply(
  seq(500, quantile(yc, .9), length.out = 8), fit_gpd_at)))
write_csv(shape_stab, file.path(VAL_DIR, "D_shape_stability.csv"))

# D3. Posterior predictive of the MAXIMUM under each body family (can the body
#     reproduce the observed extreme of 11,429?). Tail underdispersion here is a
#     red flag for extrapolating to 20,000.
pp_max <- map_dfr(c("weibull","lognormal","gamma"), function(fam) {
  fit_f <- if (fam == PRIMARY) fit else fit_hurdle(dat, fam, models$hurdle)
  post <- rstan::extract(fit_f)
  n_pos <- length(dat$pos_idx)
  vid_pos <- dat$d$vid[dat$pos_idx]
  reps <- sapply(sample(seq_along(post$b0_c), 500), function(s) {
    eta <- post$b0_c[s] + post$f_c[s, vid_pos]
    sim <- if (fam == "lognormal") rlnorm(n_pos, eta, post$sigma_ln[s])
           else if (fam == "weibull") rweibull(n_pos, post$shape[s], exp(eta))
           else rgamma(n_pos, shape = post$shape[s], rate = post$shape[s]/exp(eta))
    max(sim)
  })
  tibble(family = fam, sim_max = reps)
})
ggplot(pp_max, aes(sim_max, fill = family)) +
  geom_density(alpha = .4) +
  geom_vline(xintercept = max(yc), linetype = 2) +
  scale_x_log10() +
  labs(title = "Posterior-predictive maximum vs observed (dashed)",
       x = "simulated max (log scale)") + theme_minimal()
ggsave(file.path(VAL_DIR, "D_ppmax.png"), width = 6, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# E. Residual spatial autocorrelation (Moran's I + empirical variogram)
# -----------------------------------------------------------------------------
log_msg("E. Residual spatial autocorrelation...")
post <- rstan::extract(fit)
# village-level Pearson residual for prevalence
pi_v  <- plogis(mean(post$b0_p) + colMeans(post$f_p))
obs_v <- dat$d |> group_by(vid) |> summarise(p = mean(y_binary), n = n(), .groups="drop")
res_v <- (obs_v$p - pi_v[obs_v$vid]) * sqrt(obs_v$n)
coords_v <- as.matrix(dat$villages[, c("utmx","utmy")])

morans_I <- function(res, coords, dmin = 1) {
  # dmin (km) floors the distance so near-duplicate village sites can't blow up
  # the 1/D weights (that produced the out-of-range I ~ -35). Standard inverse-
  # distance Moran's I otherwise. S0 = sum of weights.
  D <- as.matrix(dist(coords)); D[D < dmin] <- dmin
  W <- 1 / D; diag(W) <- 0
  n <- length(res); z <- res - mean(res)
  (n / sum(W)) * (sum(W * outer(z, z)) / sum(z^2))
}
metrics$morans_I <- morans_I(res_v, coords_v)
log_msg("  Moran's I (prevalence residuals):", round(metrics$morans_I, 4))

# empirical variogram of residuals
D <- as.matrix(dist(coords_v))
vg <- map_dfr(seq(5, 120, by = 10), function(h) {
  idx <- which(D > h - 5 & D <= h + 5, arr.ind = TRUE)
  tibble(h = h, semivar = 0.5 * mean((res_v[idx[,1]] - res_v[idx[,2]])^2), n = nrow(idx))
})
ggplot(vg, aes(h, semivar)) + geom_point() + geom_line() +
  labs(title = "Residual semivariogram (flat = no leftover spatial signal)",
       x = "distance (km)", y = "semivariance") + theme_minimal()
ggsave(file.path(VAL_DIR, "E_variogram.png"), width = 6, height = 4, dpi = 150)

# -----------------------------------------------------------------------------
# F. Model comparison (LOO/WAIC) — body families (run_highrisk does the full grid)
# -----------------------------------------------------------------------------
log_msg("F. LOO/WAIC comparison...")
loo_tbl <- map_dfr(c("weibull","lognormal","gamma"), function(fam) {
  fit_f <- if (fam == PRIMARY) fit else fit_hurdle(dat, fam, models$hurdle)
  ll <- loo::extract_log_lik(fit_f, "log_lik_int", merge_chains = FALSE)
  l  <- loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
  tibble(family = fam, looic = l$estimates["looic","Estimate"],
         p_loo = l$estimates["p_loo","Estimate"],
         bad_k = mean(loo::pareto_k_values(l) > 0.7))
})
write_csv(loo_tbl, file.path(VAL_DIR, "F_loo_comparison.csv"))
print(loo_tbl)

# -----------------------------------------------------------------------------
# Save metrics summary
# -----------------------------------------------------------------------------
saveRDS(metrics, file.path(VAL_DIR, "metrics.rds"))
jsonlite::write_json(metrics, file.path(VAL_DIR, "metrics.json"),
                     auto_unbox = TRUE, pretty = TRUE)
log_msg("Validation complete. See", VAL_DIR)
