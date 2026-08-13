# =============================================================================
# NGA_shape_stability.R
# Threshold-stability scan of the GPD tail over u = 500..1500 mf/mL, so that the
# diagnostic BRACKETS the headline threshold u = 1000 (the earlier scan stopped at
# the 90th percentile, ~607, below u). Headline prior: xi ~ N(0,0.5) truncated at
# xi > -1/2, sigma ~ half-N(0,2000).
#
# Reports, per u: number of exceedances M, shape xi (median, 95% CrI), scale
# sigma_u, the threshold-stability modified scale sigma* = sigma_u - xi*u (which
# should be ~constant in u above a valid threshold), and the conditional tail
# survival P(c>20000 | c>u).
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(rstan); library(dplyr); library(ggplot2) })
VAL_DIR <- file.path(CFG$paths$out_dir, "validation")

dat <- load_nga()
yc  <- dat$d$y_intensity[dat$pos_idx]
Tt  <- CFG$hr_threshold                 # 20000

stan_gpd <- "
functions {
  real gpd_lpdf(vector e, real sigma, real xi) {
    int n = num_elements(e);
    if (fabs(xi) < 1e-8) return -n*log(sigma) - sum(e)/sigma;
    else { vector[n] z = 1 + xi*e/sigma;
      for (i in 1:n) if (z[i] <= 0) return negative_infinity();
      return -n*log(sigma) - (1+1/xi)*sum(log(z)); }
  }
}
data { int<lower=1> M; vector<lower=0>[M] e; }
parameters { real<lower=0> sigma; real<lower=-0.5> xi; }
model { sigma ~ normal(0,2000); xi ~ normal(0,0.5); e ~ gpd(sigma, xi); }
"
mod <- stan_model(model_code = stan_gpd)
qf  <- function(x) quantile(x, c(.5, .025, .975))

us  <- seq(500, 1500, by = 100)
res <- lapply(us, function(u) {
  e <- yc[yc > u] - u; M <- length(e)
  fit <- sampling(mod, data = list(M = M, e = as.array(e)),
                  chains = 4, iter = 4000, warmup = 2000, seed = CFG$seed,
                  control = list(adapt_delta = 0.999, max_treedepth = 14), refresh = 0)
  po  <- rstan::extract(fit)
  base <- 1 + po$xi*(Tt - u)/po$sigma            # GPD survival base; <=0 => beyond finite endpoint
  ps  <- ifelse(abs(po$xi) < 1e-8, exp(-(Tt - u)/po$sigma),
                ifelse(base > 0, base^(-1/po$xi), 0))
  ms  <- po$sigma - po$xi * u                      # modified (threshold-invariant) scale
  data.frame(u = u, M = M,
             xi = qf(po$xi)[1], xi_lo = qf(po$xi)[2], xi_hi = qf(po$xi)[3],
             sigma = qf(po$sigma)[1],
             modscale = qf(ms)[1], modscale_lo = qf(ms)[2], modscale_hi = qf(ms)[3],
             psurv = qf(ps)[1], psurv_lo = qf(ps)[2], psurv_hi = qf(ps)[3],
             rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE))
})
out <- do.call(rbind, res)
write.csv(out, file.path(VAL_DIR, "shape_stability_500_1500.csv"), row.names = FALSE)

# xi vs u with 95% CrI, headline u=1000 marked
p <- ggplot(out, aes(u, xi)) +
  geom_ribbon(aes(ymin = xi_lo, ymax = xi_hi), alpha = .15) +
  geom_line() + geom_point() +
  geom_vline(xintercept = 1000, linetype = 2, colour = "red") +
  labs(x = "threshold u (mf/mL)", y = expression("GPD shape "*xi),
       title = "Threshold stability of the GPD shape (500-1500 mf/mL)",
       subtitle = "dashed red = headline u = 1000") +
  theme_minimal(base_size = 11)
ggsave(file.path(VAL_DIR, "fig_shape_stability.png"), p, width = 6.5, height = 4, dpi = 150)

cat("\n=== Threshold-stability scan u = 500..1500 ===\n")
print(out, digits = 3, row.names = FALSE)
cat(sprintf("\nHeadline (u=1000): xi = %.2f [%.2f, %.2f], M = %d\n",
            out$xi[out$u==1000], out$xi_lo[out$u==1000], out$xi_hi[out$u==1000],
            out$M[out$u==1000]))
cat(sprintf("xi range across the whole scan: %.2f to %.2f\n", min(out$xi), max(out$xi)))
