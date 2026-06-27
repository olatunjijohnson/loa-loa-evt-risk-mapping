# =============================================================================
# NGA_spatial_tail_prototype.R  (Plan 2: cheap two-stage spatial-GPD prototype)
# -----------------------------------------------------------------------------
# Peter's suggestion: replace the single GPD scale sigma_u by a spatially-varying
#   sigma_u(x) = exp(alpha + beta * R(x)),   R(x) = intensity spatial field f_c(x)
# beta = 0 recovers the current pooled-GPD model.
#
# Cheap two-stage version: take the POSTERIOR MEAN of f_c at each exceedance
# village from the already-fitted lognormal hurdle (fits.rds), then fit a small
# Stan GPD with a covariate in the log-scale. Reports the posterior of beta and a
# LOO comparison of the spatial (beta free) vs pooled (beta = 0) tail.
#
# Output: results/highrisk_20000/spatial_tail/prototype_*.{csv,rds,txt}
# Runtime: a couple of minutes.
# =============================================================================
suppressMessages({
  library(rstan)
  library(loo)
})
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

ROOT <- Sys.getenv("LOA_ROOT", unset = getwd())  # run from the repository root, or set LOA_ROOT
setwd(ROOT)
source("code/NGA_common.R")

OUT <- file.path(CFG$paths$out_dir, "spatial_tail")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

u <- CFG$gpd_u   # 1000

# -----------------------------------------------------------------------------
# 1. Load the existing fit and recover f_c posterior-mean per village
# -----------------------------------------------------------------------------
log_msg("Loading fits.rds (this is large) ...")
obj  <- readRDS(file.path(CFG$paths$out_dir, "fits.rds"))
dat  <- obj$dat
post <- rstan::extract(obj$fits[["lognormal"]])   # f_c: S x V, b0_c, beta_c, sigma_ln
gpd_post <- rstan::extract(obj$gpd)               # pooled sigma, xi (current model)

# village-level posterior mean of the intensity spatial field, centred
fc_vill <- colMeans(post$f_c)            # length V
fc_vill <- fc_vill - mean(fc_vill)       # centre so alpha ~ log of pooled scale

# full intensity linear-predictor spatial part option (eta_c) — covariates at villages
# (kept for a sensitivity variant; default R = f_c)
has_X <- (!is.null(post$beta_c)) && ncol(dat$X) > 0

# -----------------------------------------------------------------------------
# 2. Exceedances over u and their village field values
# -----------------------------------------------------------------------------
yc   <- dat$d$y_intensity[dat$pos_idx]
vid  <- dat$d$vid[dat$pos_idx]
exc  <- yc > u
e    <- yc[exc] - u                      # exceedance magnitudes
Rv   <- fc_vill[vid[exc]]                # R(x) = f_c at each exceedance's village
M    <- length(e)
log_msg("Exceedances over u =", u, ":", M, "across", length(unique(vid[exc])), "villages")
log_msg("R (f_c) at exceedances: range", round(min(Rv),3), "to", round(max(Rv),3))

# -----------------------------------------------------------------------------
# 3. Stan model: spatial GPD with covariate in the log-scale
#    sigma_k = exp(alpha + beta * R_k),  e_k ~ GPD(sigma_k, xi)
# -----------------------------------------------------------------------------
STAN_SPATIAL_GPD <- "
data {
  int<lower=1> M;
  vector<lower=0>[M] e;
  vector[M] R;
  int<lower=0,upper=1> use_beta;   // 1 = spatial (beta free), 0 = pooled (beta=0)
}
parameters {
  real alpha;
  real<lower=-0.5> xi;
  array[use_beta] real beta;       // present only when use_beta = 1
}
transformed parameters {
  vector[M] sigma;
  {
    real b = use_beta ? beta[1] : 0;
    for (k in 1:M) sigma[k] = exp(alpha + b * R[k]);
  }
}
model {
  alpha ~ normal(7, 3);            // exp(7) ~ 1100, near pooled sigma_u scale
  xi ~ normal(0, 0.5);
  if (use_beta) beta[1] ~ normal(0, 1);
  for (k in 1:M) {
    real z = 1 + xi * e[k] / sigma[k];
    if (z <= 0) target += negative_infinity();
    else target += -log(sigma[k]) - (1 + 1/xi) * log(z);
  }
}
generated quantities {
  vector[M] log_lik;
  for (k in 1:M) {
    real z = 1 + xi * e[k] / sigma[k];
    log_lik[k] = (z <= 0) ? negative_infinity()
                          : -log(sigma[k]) - (1 + 1/xi) * log(z);
  }
}
"
sm <- stan_model(model_code = STAN_SPATIAL_GPD)

fit_one <- function(use_beta) {
  sampling(sm, data = list(M = M, e = e, R = Rv, use_beta = use_beta),
           chains = 4, iter = 4000, warmup = 2000, seed = CFG$seed,
           control = list(adapt_delta = 0.95), refresh = 0)
}

log_msg("Fitting POOLED tail (beta = 0) ...")
fit0 <- fit_one(0L)
log_msg("Fitting SPATIAL tail (beta free) ...")
fit1 <- fit_one(1L)

# -----------------------------------------------------------------------------
# 4. Summaries + LOO comparison
# -----------------------------------------------------------------------------
s0 <- rstan::extract(fit0); s1 <- rstan::extract(fit1)
beta_draws <- s1$beta[, 1]
beta_q <- quantile(beta_draws, c(.025,.05,.25,.5,.75,.95,.975))
p_beta_gt0 <- mean(beta_draws > 0)

ll0 <- loo::extract_log_lik(fit0, "log_lik", merge_chains = FALSE)
ll1 <- loo::extract_log_lik(fit1, "log_lik", merge_chains = FALSE)
loo0 <- loo::loo(ll0, r_eff = loo::relative_eff(exp(ll0)))
loo1 <- loo::loo(ll1, r_eff = loo::relative_eff(exp(ll1)))
cmp  <- loo::loo_compare(list(pooled = loo0, spatial = loo1))

# tail heaviness implied at the lowest vs highest f_c village (sigma ratio)
sig_lo <- exp(median(s1$alpha) + median(beta_draws) * min(Rv))
sig_hi <- exp(median(s1$alpha) + median(beta_draws) * max(Rv))

sink(file.path(OUT, "prototype_summary.txt"))
cat("Spatial-GPD prototype (two-stage, plug-in posterior-mean f_c)\n")
cat("=============================================================\n\n")
cat(sprintf("Exceedances over u=%d: M=%d across %d villages\n\n",
            u, M, length(unique(vid[exc]))))
cat("POOLED model (beta=0, ~ current model):\n")
cat(sprintf("  alpha (log sigma_u)  median %.3f  -> sigma_u ~ %.0f\n",
            median(s0$alpha), exp(median(s0$alpha))))
cat(sprintf("  xi                   median %.3f  [%.3f, %.3f]\n\n",
            median(s0$xi), quantile(s0$xi,.025), quantile(s0$xi,.975)))
cat("SPATIAL model (beta free):\n")
cat(sprintf("  alpha   median %.3f\n", median(s1$alpha)))
cat(sprintf("  xi      median %.3f  [%.3f, %.3f]\n", median(s1$xi),
            quantile(s1$xi,.025), quantile(s1$xi,.975)))
cat(sprintf("  beta    median %.3f  90%% [%.3f, %.3f]  95%% [%.3f, %.3f]\n",
            beta_q["50%"], beta_q["5%"], beta_q["95%"], beta_q["2.5%"], beta_q["97.5%"]))
cat(sprintf("  P(beta > 0) = %.3f\n", p_beta_gt0))
cat(sprintf("  implied GPD scale: low-intensity village %.0f  vs  high-intensity village %.0f  (ratio %.2fx)\n\n",
            sig_lo, sig_hi, sig_hi/sig_lo))
cat("LOO comparison (does beta earn its keep?):\n")
print(cmp)
cat("\nInterpretation: if elpd_diff for 'spatial' is not clearly larger than its se_diff,\n")
cat("the extra parameter does NOT earn its keep on this dataset (as Peter anticipated).\n")
sink()

cat(readLines(file.path(OUT, "prototype_summary.txt")), sep = "\n")

saveRDS(list(fit0 = fit0, fit1 = fit1, beta = beta_draws, cmp = cmp,
             M = M, Rv = Rv, e = e),
        file.path(OUT, "prototype_fits.rds"))
write.csv(data.frame(beta = beta_draws), file.path(OUT, "prototype_beta_draws.csv"),
          row.names = FALSE)
log_msg("Prototype done. Output in", OUT)
