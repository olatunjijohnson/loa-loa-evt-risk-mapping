# =============================================================================
# NGA_prior_sensitivity_xi.R
# Reviewer C: prior sensitivity of the GPD tail shape xi.
# With only M=37 exceedances over u=1000 mf/mL the shape is informed partly by
# the prior, so we refit the pooled GPD tail under several priors on xi and show
# that the quantity that actually matters --- the conditional tail survival
# P(c > 20000 | c > u) --- is stable. Because the per-LGA high-risk burden is
# (body term P(c>u|x)) x (this conditional survival), stability here implies the
# per-LGA conclusion is insensitive to the xi prior.
#
# Outputs:
#   results/highrisk_20000/validation/prior_sensitivity_xi.csv
#   results/highrisk_20000/validation/fig_prior_xi.png
#     (also mirrored to $LOA_MS_FIGDIR for the SI, if that variable is set)
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

VAL_DIR <- file.path(CFG$paths$out_dir, "validation")
dir.create(VAL_DIR, recursive = TRUE, showWarnings = FALSE)

dat <- load_nga()
yc  <- dat$d$y_intensity[dat$pos_idx]
u   <- CFG$gpd_u                       # 1000
e   <- yc[yc > u] - u
M   <- length(e)
Tt  <- CFG$hr_threshold                # 20000
cat(sprintf("Exceedances over u=%d: M=%d (max density %.0f)\n", u, M, max(yc)))

# Parametrised GPD: prior sd on xi passed as data; lower truncation fixed at
# compile time, so we build one model per truncation bound.
stan_gpd <- function(xi_lb) sprintf("
functions {
  real gpd_lpdf(vector e, real sigma, real xi) {
    int n = num_elements(e);
    if (fabs(xi) < 1e-8) return -n*log(sigma) - sum(e)/sigma;
    else {
      vector[n] z = 1 + xi * e / sigma;
      for (i in 1:n) if (z[i] <= 0) return negative_infinity();
      return -n*log(sigma) - (1 + 1/xi) * sum(log(z));
    }
  }
}
data { int<lower=1> M; vector<lower=0>[M] e; real<lower=0> prior_xi_sd; }
parameters { real<lower=0> sigma; real<lower=%s> xi; }
model { sigma ~ normal(0, 2000); xi ~ normal(0, prior_xi_sd); e ~ gpd(sigma, xi); }
", format(xi_lb, nsmall = 1))

m_lb05 <- stan_model(model_code = stan_gpd(-0.5))
m_lb10 <- stan_model(model_code = stan_gpd(-1.0))

# Conditional tail survival P(c>T | c>u) per posterior draw.
cond_surv <- function(xi, sigma, Tu) {
  ifelse(abs(xi) < 1e-8, exp(-Tu / sigma), (1 + xi * Tu / sigma)^(-1 / xi))
}

scenarios <- tribble(
  ~label,                 ~prior_xi_sd, ~model,
  "N(0,0.5), xi>-0.5 (headline)", 0.5,  "lb05",
  "N(0,1),   xi>-0.5",            1.0,  "lb05",
  "N(0,2),   xi>-0.5",            2.0,  "lb05",
  "N(0,0.5), xi>-1.0",           0.5,  "lb10"
)

res <- list(); draws <- list()
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  mod <- if (s$model == "lb05") m_lb05 else m_lb10
  fit <- sampling(mod, data = list(M = M, e = as.array(e), prior_xi_sd = s$prior_xi_sd),
                  chains = 4, iter = 4000, warmup = 2000, seed = CFG$seed,
                  control = list(adapt_delta = 0.999, max_treedepth = 14), refresh = 0)
  po <- rstan::extract(fit)
  ps <- cond_surv(po$xi, po$sigma, Tt - u)
  q  <- function(x) quantile(x, c(.5, .025, .975))
  res[[i]] <- tibble(
    scenario = s$label, M = M,
    xi_med = q(po$xi)[1],      xi_lo = q(po$xi)[2],      xi_hi = q(po$xi)[3],
    sigma_med = q(po$sigma)[1], sigma_lo = q(po$sigma)[2], sigma_hi = q(po$sigma)[3],
    psurv_med = q(ps)[1],      psurv_lo = q(ps)[2],      psurv_hi = q(ps)[3],
    max_rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE)
  )
  draws[[i]] <- tibble(scenario = s$label, xi = po$xi, psurv = ps)
  cat(sprintf("  %-28s xi=%.2f [%.2f,%.2f]  P(c>20000|c>u)=%.2e [%.2e,%.2e]\n",
              s$label, res[[i]]$xi_med, res[[i]]$xi_lo, res[[i]]$xi_hi,
              res[[i]]$psurv_med, res[[i]]$psurv_lo, res[[i]]$psurv_hi))
}
out <- bind_rows(res)
write.csv(out, file.path(VAL_DIR, "prior_sensitivity_xi.csv"), row.names = FALSE)

# Figure: posterior of the conditional tail survival under each prior (log10).
dd <- bind_rows(draws) |> mutate(scenario = factor(scenario, levels = scenarios$label))
p <- ggplot(dd, aes(x = pmax(psurv, 1e-12), colour = scenario)) +
  geom_density(linewidth = 0.7) +
  scale_x_log10() +
  labs(x = expression(P(c > 20000 ~ "|" ~ c > u)),
       y = "posterior density", colour = NULL,
       title = "Prior sensitivity of the GPD tail (M = 37 exceedances)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(VAL_DIR, "fig_prior_xi.png"), p, width = 7, height = 4.2, dpi = 150)
# Optionally mirror the figure into a manuscript figures directory:
#   LOA_MS_FIGDIR=/path/to/manuscript/figures Rscript code/NGA_prior_sensitivity_xi.R
ms_figdir <- Sys.getenv("LOA_MS_FIGDIR", unset = "")
if (nzchar(ms_figdir) && dir.exists(ms_figdir))
  file.copy(file.path(VAL_DIR, "fig_prior_xi.png"),
            file.path(ms_figdir, "fig_prior_xi.png"), overwrite = TRUE)
cat("\nWrote prior_sensitivity_xi.csv and fig_prior_xi.png\n")
print(out)
