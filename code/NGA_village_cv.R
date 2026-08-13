# =============================================================================
# NGA_village_cv.R
# Reviewer D: village-level (leave-one-village-out style) cross-validation.
# The block CV in NGA_validation.R removes whole spatial regions (k-means folds)
# and is deliberately harsh (held-out AUC ~ 0.52). Here we use RANDOM folds, so
# held-out villages are interspersed among retained neighbours: every village is
# predicted from a model fitted WITHOUT it, but local neighbours remain. This
# measures LOCAL interpolation skill at typical inter-village distances.
#
# Same model, same kriging predictor, same metrics as the block CV, so the two
# are directly comparable.
#
# Outputs:
#   results/highrisk_20000/validation/village_cv_metrics.json
#   results/highrisk_20000/validation/fig_village_cv.png
#     (also mirrored to $LOA_MS_FIGDIR for the SI, if that variable is set)
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(loo); library(pROC); library(scoringRules)
                   library(ggplot2); library(dplyr); library(purrr); library(jsonlite) })

VAL_DIR <- file.path(CFG$paths$out_dir, "validation")
dir.create(VAL_DIR, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

PRIMARY <- "lognormal"
K       <- as.integer(Sys.getenv("VILLAGE_CV_K", "10"))   # random folds
CV_ITER <- as.integer(Sys.getenv("VILLAGE_CV_ITER", "800"))
CV_WARM <- as.integer(Sys.getenv("VILLAGE_CV_WARM", "400"))

dat    <- load_nga()
models <- compile_models()

# RANDOM folds (interspersed), in contrast to the k-means spatial blocks.
make_folds_random <- function(villages, K = 10, seed = 1) {
  set.seed(seed)
  villages$fold <- sample(rep_len(seq_len(K), nrow(villages)))
  villages
}

# Kriging predictor for held-out villages (verbatim from NGA_validation.R).
predict_holdout <- function(fit_tr, dat_tr, hold_coords, X_hold, family) {
  post  <- rstan::extract(fit_tr)
  Dgv_h <- fields_rdist(hold_coords, as.matrix(dat_tr$villages[, c("utmx", "utmy")]))
  f_p_h <- krige_mean(post$f_p, dat_tr$Dmat, Dgv_h, post$rho)
  f_c_h <- krige_mean(post$f_c, dat_tr$Dmat, Dgv_h, post$rho)
  eta_p_h <- sweep(f_p_h, 1, post$b0_p, "+")
  eta_c_h <- sweep(f_c_h, 1, post$b0_c, "+")
  if (ncol(X_hold) > 0) {
    eta_p_h <- eta_p_h + post$beta_p %*% t(X_hold)
    eta_c_h <- eta_c_h + post$beta_c %*% t(X_hold)
  }
  list(pi = plogis(eta_p_h), eta_c = eta_c_h,
       shape = post$shape, sigma_ln = post$sigma_ln)
}

# CV driver (verbatim logic from run_spatial_cv, random folds).
run_cv <- function(dat, family, folds, iter, warmup) {
  preds <- list()
  for (k in sort(unique(folds$fold))) {
    log_msg("  fold", k, "/", max(folds$fold))
    tr_vill <- folds |> filter(fold != k)
    dtr <- dat$d |> filter(vill_key %in% tr_vill$vill_key)
    tr_vill <- tr_vill |> mutate(vid = row_number())
    dtr <- dtr |> select(-vid) |> left_join(tr_vill |> select(vill_key, vid), by = "vill_key")
    Xtr <- if (ncol(dat$Xv) > 0) dat$Xv[dtr$vill_key, , drop = FALSE] else matrix(0, nrow(dtr), 0)
    dat_tr <- list(d = dtr, villages = tr_vill, V = nrow(tr_vill),
                   Dmat = as.matrix(dist(as.matrix(tr_vill[, c("utmx","utmy")]))),
                   X = Xtr, pos_idx = which(dtr$y_binary == 1))
    fit_tr <- fit_hurdle(dat_tr, family, models$hurdle, iter = iter, warmup = warmup)
    ho_vill <- folds |> filter(fold == k)
    ho_ind  <- dat$d |> filter(vill_key %in% ho_vill$vill_key)
    Xho <- if (ncol(dat$Xv) > 0) dat$Xv[ho_vill$vill_key, , drop = FALSE] else matrix(0, nrow(ho_vill), 0)
    p <- predict_holdout(fit_tr, dat_tr,
                         as.matrix(ho_vill[, c("utmx","utmy")]), Xho, family)
    vmap <- match(ho_ind$vill_key, ho_vill$vill_key)
    preds[[as.character(k)]] <- list(ind = ho_ind, vmap = vmap, p = p)
  }
  preds
}

log_msg("Village-level random", K, "-fold CV (", CV_ITER, "iter ) ...")
folds <- make_folds_random(dat$villages, K, seed = CFG$seed)
cv <- run_cv(dat, PRIMARY, folds, CV_ITER, CV_WARM)

# --- Prevalence metrics ---
yb     <- unlist(lapply(cv, function(f) f$ind$y_binary))
pi_hat <- unlist(lapply(cv, function(f) colMeans(f$p$pi)[f$vmap]))
roc_obj <- pROC::roc(yb, pi_hat, quiet = TRUE)
prev <- list(auc = as.numeric(pROC::auc(roc_obj)),
             brier = mean((pi_hat - yb)^2),
             logloss = -mean(yb*log(pi_hat) + (1-yb)*log(1-pi_hat)))

# --- Intensity metrics (CRPS, 90% coverage) ---
intensity <- map_dfr(cv, function(f) {
  pos <- which(f$ind$y_binary == 1)
  if (!length(pos)) return(tibble())
  map_dfr(pos, function(j) {
    v <- f$vmap[j]; eta <- f$p$eta_c[, v]
    draws <- rlnorm(length(eta), eta, f$p$sigma_ln)
    obs <- f$ind$y_intensity[j]
    tibble(crps = scoringRules::crps_sample(obs, draws),
           in90 = obs >= quantile(draws, .05) & obs <= quantile(draws, .95))
  })
})
out <- list(method = "random village-level k-fold CV", K = K,
            iter = CV_ITER, warmup = CV_WARM,
            n_villages = dat$V, n_individuals = nrow(dat$d),
            prevalence = prev,
            intensity = list(mean_crps = mean(intensity$crps),
                             cover90 = mean(intensity$in90)),
            block_cv_reference = list(auc = 0.5225, cover90 = 0.922))
write_json(out, file.path(VAL_DIR, "village_cv_metrics.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 4)

# --- Reliability figure (random village CV vs the harsh block CV) ---
rel <- tibble(p = pi_hat, y = yb) |>
  mutate(bin = cut(p, seq(0, 1, 0.1), include.lowest = TRUE)) |>
  group_by(bin) |> summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop")
g <- ggplot(rel, aes(pred, obs, size = n)) +
  geom_abline(linetype = 2) + geom_point(alpha = .8) +
  labs(title = sprintf("Prevalence reliability, random village %d-fold CV (AUC %.2f)",
                       K, prev$auc),
       x = "Predicted prevalence", y = "Observed prevalence") +
  theme_minimal(base_size = 11)
ggsave(file.path(VAL_DIR, "fig_village_cv.png"), g, width = 5, height = 5, dpi = 150)
# Optionally mirror the figure into a manuscript figures directory:
#   LOA_MS_FIGDIR=/path/to/manuscript/figures Rscript code/NGA_village_cv.R
ms_figdir <- Sys.getenv("LOA_MS_FIGDIR", unset = "")
if (nzchar(ms_figdir) && dir.exists(ms_figdir))
  file.copy(file.path(VAL_DIR, "fig_village_cv.png"),
            file.path(ms_figdir, "fig_village_cv.png"), overwrite = TRUE)

cat(sprintf("\n=== Village-level random %d-fold CV ===\n", K))
cat(sprintf("Prevalence : AUC=%.3f  Brier=%.4f  logloss=%.4f\n",
            prev$auc, prev$brier, prev$logloss))
cat(sprintf("Intensity  : mean CRPS=%.1f  90%% coverage=%.3f\n",
            mean(intensity$crps), mean(intensity$in90)))
cat(sprintf("(block CV reference: AUC=0.52, coverage=0.92)\n"))
