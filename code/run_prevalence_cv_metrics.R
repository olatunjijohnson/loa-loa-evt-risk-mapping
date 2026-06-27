# =============================================================================
# run_prevalence_cv_metrics.R — extra prevalence-component metrics from the
# spatial block CV: PR-AUC (vs base-rate baseline) and Brier / log-loss skill
# scores. Re-runs only the 6-fold spatial CV (no full-validation rerun).
# Functions copied verbatim from NGA_validation.R so the CV is identical.
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(dplyr); library(pROC) })
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")
PRIMARY <- "lognormal"

dat <- load_nga(); models <- compile_models()

make_folds <- function(villages, K = 10, seed = 1) {
  set.seed(seed)
  km <- kmeans(as.matrix(villages[, c("utmx", "utmy")]), centers = K, nstart = 5)
  villages$fold <- km$cluster; villages
}
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
  list(pi = plogis(eta_p_h), eta_c = eta_c_h, shape = post$shape, sigma_ln = post$sigma_ln)
}
run_spatial_cv <- function(dat, family, K = 6, iter = 800, warmup = 400) {
  vill <- make_folds(dat$villages, K); preds <- list()
  for (k in seq_len(K)) {
    log_msg("  CV fold", k, "/", K)
    tr_vill <- vill |> filter(fold != k)
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
    p <- predict_holdout(fit_tr, dat_tr, as.matrix(ho_vill[, c("utmx","utmy")]), Xho, family)
    vmap <- match(ho_ind$vill_key, ho_vill$vill_key)
    preds[[k]] <- list(ind = ho_ind, vmap = vmap, p = p)
  }
  preds
}

log_msg("Spatial block CV (6 folds)...")
cv <- run_spatial_cv(dat, PRIMARY, K = 6)
yb <- unlist(lapply(cv, function(f) f$ind$y_binary))
pi_hat <- unlist(lapply(cv, function(f) colMeans(f$p$pi)[f$vmap]))
pc <- pmin(pmax(pi_hat, 1e-6), 1 - 1e-6)             # clamp for log-loss
base <- mean(yb)

roc_auc <- as.numeric(pROC::auc(pROC::roc(yb, pi_hat, quiet = TRUE)))
pr_auc <- if (requireNamespace("PRROC", quietly = TRUE)) {
  PRROC::pr.curve(scores.class0 = pi_hat[yb == 1], scores.class1 = pi_hat[yb == 0])$auc.integral
} else {                                              # average precision (manual)
  o <- order(pi_hat, decreasing = TRUE); y <- yb[o]
  tp <- cumsum(y); prec <- tp / seq_along(y); rec <- tp / sum(y)
  sum(prec * diff(c(0, rec)))
}
brier_m <- mean((pi_hat - yb)^2); brier_b <- mean((base - yb)^2); bss <- 1 - brier_m / brier_b
ll_m <- -mean(yb * log(pc) + (1 - yb) * log(1 - pc))
ll_b <- -mean(yb * log(base) + (1 - yb) * log(1 - base)); llss <- 1 - ll_m / ll_b

cat("\n================  PREVALENCE SPATIAL BLOCK-CV  ================\n")
cat(sprintf("n individuals        : %d   (positives %d, base rate %.4f)\n", length(yb), sum(yb), base))
cat(sprintf("ROC-AUC              : %.3f   (no-skill 0.500)\n", roc_auc))
cat(sprintf("PR-AUC (avg prec.)   : %.3f   (no-skill = base rate %.3f)\n", pr_auc, base))
cat(sprintf("Brier  model/base    : %.4f / %.4f   -> BSS  %.3f\n", brier_m, brier_b, bss))
cat(sprintf("LogLoss model/base   : %.4f / %.4f   -> skill %.3f\n", ll_m, ll_b, llss))
cat("==============================================================\n")
jsonlite::write_json(list(n = length(yb), positives = sum(yb), base_rate = base,
  roc_auc = roc_auc, pr_auc = pr_auc, brier_model = brier_m, brier_baseline = brier_b,
  bss = bss, logloss_model = ll_m, logloss_baseline = ll_b, logloss_skill = llss),
  "results/highrisk_20000/validation/prevalence_cv_extra.json", auto_unbox = TRUE, pretty = TRUE)
