# =============================================================================
# NGA_joint_tail.R   (Plan 3: full joint GPD-in-hurdle model)
# -----------------------------------------------------------------------------
# Peter's reformulation, fitted JOINTLY (uncertainty in f_c propagated):
#   lognormal hurdle body  +  spatial GPD tail with
#       sigma_u(x) = exp(alpha + beta * f_c(x)),   xi shared
#   beta = 0  -> the current pooled-GPD model (structurally)
#   beta free -> tail heaviness tracks the intensity spatial field
#
# Fits BOTH variants (est_beta = 1 and 0), recomputes the per-LGA answer
# P(q_g < 0.1%) under the spatially-varying tail, and writes a LOO comparison
# (does beta earn its keep?) plus a side-by-side vs the committed headline.
#
# Env vars:
#   JOINT_SMOKE=1  -> tiny fit (2 chains, 200 iter) to check it compiles/runs.
# Output: results/highrisk_20000/spatial_tail/joint_*.{rds,csv,txt}
# =============================================================================
suppressMessages({ library(rstan); library(loo); library(dplyr); library(purrr); library(readr) })
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ROOT <- Sys.getenv("LOA_ROOT", unset = getwd())  # run from the repository root, or set LOA_ROOT
setwd(ROOT)
source("code/NGA_common.R")

SMOKE <- nzchar(Sys.getenv("JOINT_SMOKE"))
OUT   <- file.path(CFG$paths$out_dir, "spatial_tail")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")
u <- CFG$gpd_u

# -----------------------------------------------------------------------------
# Joint Stan model: lognormal hurdle body + spatial GPD tail.
# (Body block identical to STAN_HURDLE; only the tail data/params/likelihood added.)
# -----------------------------------------------------------------------------
STAN_HURDLE_TAIL <- "
data {
  int<lower=1> N;
  int<lower=1> V;
  array[N] int<lower=1,upper=V> vid;
  array[N] int<lower=0,upper=1> y;
  int<lower=0> Npos;
  array[Npos] int<lower=1,upper=N> pos_idx;
  vector<lower=0>[Npos] yc;
  int<lower=0> P;
  matrix[N, P] X;
  matrix[V, V] Dmat;
  real<lower=0> hs_scale_global;
  // --- spatial-tail block ---
  int<lower=1> M;                       // # exceedances over u
  array[M] int<lower=1,upper=V> exc_vid; // village of each exceedance
  vector<lower=0>[M] e;                 // exceedance magnitude (count - u)
  int<lower=0,upper=1> est_beta;        // 1 = spatial tail, 0 = pooled tail
}
transformed data {
  real jitter = 1e-8;
  real hs_slab_scale = 2;
  real hs_slab_df = 4;
}
parameters {
  real b0_p;
  real b0_c;
  vector[P] zb_p;  vector<lower=0>[P] hs_lam_p; real<lower=0> hs_tau_p; real<lower=0> hs_caux_p;
  vector[P] zb_c;  vector<lower=0>[P] hs_lam_c; real<lower=0> hs_tau_c; real<lower=0> hs_caux_c;
  real<lower=0> rho;
  vector<lower=0>[2] sigma_f;
  cholesky_factor_corr[2] L_Omega;
  matrix[V, 2] Z;
  real<lower=0> sigma_ln;
  // --- spatial-tail params ---
  real alpha_tail;                      // baseline log GPD scale
  real<lower=-0.5> xi_tail;             // GPD shape
  array[est_beta] real beta_tail;       // field->scale slope (only when est_beta=1)
}
transformed parameters {
  vector[V] f_p; vector[V] f_c; vector[P] beta_p; vector[P] beta_c;
  {
    real c2p = square(hs_slab_scale) * hs_caux_p;
    real c2c = square(hs_slab_scale) * hs_caux_c;
    vector[P] lt_p = sqrt(c2p * square(hs_lam_p) ./ (c2p + square(hs_tau_p) * square(hs_lam_p)));
    vector[P] lt_c = sqrt(c2c * square(hs_lam_c) ./ (c2c + square(hs_tau_c) * square(hs_lam_c)));
    beta_p = zb_p .* lt_p * hs_tau_p;
    beta_c = zb_c .* lt_c * hs_tau_c;
  }
  {
    matrix[V, V] Lk = cholesky_decompose(add_diag(exp(-Dmat / rho), jitter));
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_f, L_Omega);
    matrix[V, 2] F = Lk * Z * L_Sigma';
    f_p = F[, 1];  f_c = F[, 2];
  }
}
model {
  b0_p ~ normal(0, 2);  b0_c ~ normal(0, 2);
  zb_p ~ std_normal();  zb_c ~ std_normal();
  hs_lam_p ~ cauchy(0, 1);  hs_lam_c ~ cauchy(0, 1);
  hs_tau_p ~ cauchy(0, hs_scale_global);  hs_tau_c ~ cauchy(0, hs_scale_global);
  hs_caux_p ~ inv_gamma(0.5 * hs_slab_df, 0.5 * hs_slab_df);
  hs_caux_c ~ inv_gamma(0.5 * hs_slab_df, 0.5 * hs_slab_df);
  rho ~ inv_gamma(5, 150);
  sigma_f ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(Z) ~ std_normal();
  sigma_ln ~ normal(0, 1);
  // tail priors
  alpha_tail ~ normal(7, 3);            // exp(7) ~ 1100, near pooled sigma_u
  xi_tail ~ normal(0, 0.5);
  if (est_beta) beta_tail[1] ~ normal(0, 1);

  // prevalence
  { vector[N] eta_p = b0_p + f_p[vid]; if (P > 0) eta_p += X * beta_p; y ~ bernoulli_logit(eta_p); }
  // intensity body (lognormal, positives)
  for (k in 1:Npos) {
    int i = pos_idx[k];
    real eta = b0_c + f_c[vid[i]]; if (P > 0) eta += X[i] * beta_c;
    yc[k] ~ lognormal(eta, sigma_ln);
  }
  // spatial GPD tail on exceedances
  {
    real b = est_beta ? beta_tail[1] : 0;
    for (m in 1:M) {
      real sigma = exp(alpha_tail + b * f_c[exc_vid[m]]);
      real z = 1 + xi_tail * e[m] / sigma;
      if (z <= 0) target += negative_infinity();
      else target += -log(sigma) - (1 + 1/xi_tail) * log(z);
    }
  }
}
generated quantities {
  real rho_pc = multiply_lower_tri_self_transpose(L_Omega)[2, 1];
  vector[Npos] log_lik_int;
  vector[M] log_lik_tail;
  for (k in 1:Npos) {
    int i = pos_idx[k];
    real eta = b0_c + f_c[vid[i]]; if (P > 0) eta += X[i] * beta_c;
    log_lik_int[k] = lognormal_lpdf(yc[k] | eta, sigma_ln);
  }
  {
    real b = est_beta ? beta_tail[1] : 0;
    for (m in 1:M) {
      real sigma = exp(alpha_tail + b * f_c[exc_vid[m]]);
      real z = 1 + xi_tail * e[m] / sigma;
      log_lik_tail[m] = (z <= 0) ? negative_infinity()
                                 : -log(sigma) - (1 + 1/xi_tail) * log(z);
    }
  }
}
"

# -----------------------------------------------------------------------------
# Data + exceedances
# -----------------------------------------------------------------------------
log_msg("Loading data ...")
dat <- load_nga()
yc      <- dat$d$y_intensity[dat$pos_idx]
vid_pos <- dat$d$vid[dat$pos_idx]
exc     <- yc > u
e       <- yc[exc] - u
exc_vid <- vid_pos[exc]
M       <- length(e)
log_msg("N =", nrow(dat$d), "| V =", dat$V, "| positives =", length(dat$pos_idx),
        "| exceedances over u=", u, ":", M)

P  <- ncol(dat$X); Nn <- nrow(dat$d); p0 <- CFG$hs_p0
tau0 <- if (P > p0) p0 / (P - p0) / sqrt(Nn) else 1
base_data <- list(
  N = Nn, V = dat$V, vid = dat$d$vid, y = dat$d$y_binary,
  Npos = length(dat$pos_idx), pos_idx = dat$pos_idx, yc = yc,
  P = P, X = dat$X, Dmat = dat$Dmat, hs_scale_global = tau0,
  M = M, exc_vid = exc_vid, e = e)

log_msg("Compiling joint hurdle+tail model ...")
sm <- stan_model(model_code = STAN_HURDLE_TAIL)

iter   <- if (SMOKE) 200  else CFG$iter
warm   <- if (SMOKE) 100  else CFG$warmup
chains <- if (SMOKE) 2    else CFG$chains
fit_variant <- function(est_beta) {
  sampling(sm, data = c(base_data, list(est_beta = est_beta)),
           chains = chains, iter = iter, warmup = warm, seed = CFG$seed,
           control = list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth),
           refresh = if (SMOKE) 50 else 200)
}

log_msg("Fitting SPATIAL tail variant (est_beta = 1) ...")
fit_sp <- fit_variant(1L)
log_msg("Fitting POOLED tail variant (est_beta = 0) ...")
fit_po <- fit_variant(0L)
saveRDS(list(fit_sp = fit_sp, fit_po = fit_po, dat = dat,
             exc_vid = exc_vid, e = e, u = u),
        file.path(OUT, if (SMOKE) "joint_fits_smoke.rds" else "joint_fits.rds"))

# -----------------------------------------------------------------------------
# Per-cell GPD ratio with spatially-varying scale, and per-LGA q_g
# -----------------------------------------------------------------------------
gpd_ratio_vec <- function(sigma, xi, u, Tval) {
  z <- 1 + xi * (Tval - u) / sigma
  ifelse(z <= 0, 0, z^(-1 / xi))
}

# Variant of predict_lga_exact: gpd tail with sigma_u(x) = exp(alpha + beta*f_c(x)).
predict_lga_spatialtail <- function(post, dat, geom, Tval, est_beta,
                                    cfg = CFG, thin = cfg$thin_pred, seed = cfg$seed) {
  set.seed(seed)
  S0 <- length(post$b0_p); idx <- seq(1, S0, by = thin); S <- length(idx)
  Dmat <- dat$Dmat; Pc <- ncol(dat$X_grid)
  qg <- matrix(0, S, length(geom))
  for (ii in seq_along(idx)) {
    s   <- idx[ii]; rho <- post$rho[s]
    Kvv <- exp(-Dmat / rho); diag(Kvv) <- diag(Kvv) + 1e-8; Kvv_i <- solve(Kvv)
    L_Sigma <- diag(c(post$sigma_f[s, 1], post$sigma_f[s, 2])) %*% post$L_Omega[s, , ]
    b   <- if (est_beta) post$beta_tail[s] else 0
    for (g in seq_along(geom)) {
      G    <- geom[[g]]
      cd   <- .cond_draw_icm(post$f_p[s, ], post$f_c[s, ], Kvv_i, G$Dgv, G$Dgg, rho, L_Sigma)
      eta_p <- post$b0_p[s] + cd$fp
      eta_c <- post$b0_c[s] + cd$fc
      if (Pc > 0) {
        Xg <- dat$X_grid[G$cells, , drop = FALSE]
        eta_p <- eta_p + as.numeric(Xg %*% post$beta_p[s, ])
        eta_c <- eta_c + as.numeric(Xg %*% post$beta_c[s, ])
      }
      pi_g  <- plogis(eta_p)
      su    <- surv_cells("lognormal", eta_c, sigma_ln = post$sigma_ln[s], Tval = cfg$gpd_u)
      sigma_cell <- exp(post$alpha_tail[s] + b * cd$fc)        # R(x) = f_c(x)
      R     <- gpd_ratio_vec(sigma_cell, post$xi_tail[s], cfg$gpd_u, Tval)
      ix    <- pi_g * su * R
      qg[ii, g] <- weighted.mean(ix, G$pop)
    }
  }
  qg
}

log_msg("Precomputing LGA geometry ...")
geom <- precompute_lga_geometry(dat)
post_sp <- rstan::extract(fit_sp)
post_po <- rstan::extract(fit_po)

log_msg("Predicting per-LGA q_g (spatial tail) ...")
qg_sp <- predict_lga_spatialtail(post_sp, dat, geom, CFG$hr_threshold, est_beta = 1L)
log_msg("Predicting per-LGA q_g (pooled tail, joint) ...")
qg_po <- predict_lga_spatialtail(post_po, dat, geom, CFG$hr_threshold, est_beta = 0L)

sum_sp <- summarise_qg(qg_sp, geom) |> mutate(variant = "spatial_tail")
sum_po <- summarise_qg(qg_po, geom) |> mutate(variant = "pooled_tail_joint")
write_csv(bind_rows(sum_sp, sum_po), file.path(OUT, "joint_qg_summaries.csv"))

# -----------------------------------------------------------------------------
# beta posterior + LOO comparison (does beta earn its keep, jointly?)
# -----------------------------------------------------------------------------
beta_draws <- post_sp$beta_tail
ll_sp <- loo::extract_log_lik(fit_sp, "log_lik_tail", merge_chains = FALSE)
ll_po <- loo::extract_log_lik(fit_po, "log_lik_tail", merge_chains = FALSE)
loo_sp <- loo::loo(ll_sp, r_eff = loo::relative_eff(exp(ll_sp)))
loo_po <- loo::loo(ll_po, r_eff = loo::relative_eff(exp(ll_po)))
cmp <- loo::loo_compare(list(pooled = loo_po, spatial = loo_sp))

sink(file.path(OUT, if (SMOKE) "joint_summary_smoke.txt" else "joint_summary.txt"))
cat("Joint hurdle + spatial-GPD tail  (", if (SMOKE) "SMOKE" else "FULL", ")\n")
cat("=================================================\n\n")
cat(sprintf("Exceedances M = %d over u = %d\n\n", M, u))
rh <- function(fit, par) {
  sm <- rstan::summary(fit, pars = par)$summary; sm[, c("mean","2.5%","50%","97.5%","n_eff","Rhat")]
}
cat("SPATIAL variant tail params:\n"); print(round(rh(fit_sp, c("alpha_tail","beta_tail","xi_tail")), 3))
cat(sprintf("\n  P(beta > 0) = %.3f\n", mean(beta_draws > 0)))
cat("\nPOOLED variant tail params:\n"); print(round(rh(fit_po, c("alpha_tail","xi_tail")), 3))
cat("\nLOO comparison (tail log-lik) — does beta earn its keep?\n"); print(cmp)
cat("\nMax Rhat (spatial fit):", round(max(rstan::summary(fit_sp)$summary[,"Rhat"], na.rm=TRUE), 3), "\n")
cat("Per-LGA answer P(q_g<0.1%):\n")
cat(sprintf("  spatial tail : min %.4f  median %.4f  max %.4f\n",
            min(sum_sp$p_below_0.1pct), median(sum_sp$p_below_0.1pct), max(sum_sp$p_below_0.1pct)))
cat(sprintf("  pooled  (joint): min %.4f  median %.4f  max %.4f\n",
            min(sum_po$p_below_0.1pct), median(sum_po$p_below_0.1pct), max(sum_po$p_below_0.1pct)))
cat(sprintf("  worst-LGA burden per 100k: spatial %.3f  vs pooled %.3f\n",
            max(sum_sp$q_mean)*1e5, max(sum_po$q_mean)*1e5))
sink()
cat(readLines(file.path(OUT, if (SMOKE) "joint_summary_smoke.txt" else "joint_summary.txt")), sep = "\n")
log_msg("joint_tail done. Output in", OUT)
