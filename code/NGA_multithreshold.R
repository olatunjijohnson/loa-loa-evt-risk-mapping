# =============================================================================
# NGA_multithreshold.R — exceedance probabilities at 10,000 / 15,000 / 20,000 mf/mL.
# Same headline model (lognormal body + pooled GPD tail). For the GPD tail the
# threshold T enters ONLY through the per-draw ratio R(T) = gpd_ratio(.,.,u,T); the
# expensive per-LGA conditional draws (pi_g * S_lognormal(u)) are identical across
# thresholds. So we run the per-draw machinery ONCE (base = q_g with R=1) and rescale
# to each T. The T=20,000 column reproduces results/.../q_g_summaries.csv as a check.
#   Run AFTER NGA_run_highrisk.R.  Usage:  Rscript code/NGA_multithreshold.R
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(dplyr); library(purrr); library(readr); library(tidyr) })
OUT <- CFG$paths$out_dir
THR <- c(10000, 15000, 20000)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

log_msg("Loading fits + data...")
fr  <- readRDS(file.path(OUT, "fits.rds")); dat <- fr$dat
post     <- rstan::extract(fr$fits$lognormal)
gpd_post <- rstan::extract(fr$gpd)
geom <- precompute_lga_geometry(dat)

# ---- base draws (q_g with R = 1), reproducing predict_lga_exact's RNG stream ----
# Mirrors predict_lga_exact(tail = "gpd") but factors R(T) out of the cell sum so
# one pass serves all thresholds. R is a per-draw scalar (constant across an LGA's
# cells), so weighted.mean(pi_g*su*R) = R * weighted.mean(pi_g*su).
predict_lga_base <- function(post, dat, geom, cfg = CFG,
                             thin = cfg$thin_pred, seed = cfg$seed) {
  set.seed(seed)
  S0 <- length(post$b0_p); idx <- seq(1, S0, by = thin); S <- length(idx)
  Dmat <- dat$Dmat; Pc <- ncol(dat$X_grid)
  base <- matrix(0, S, length(geom))
  for (ii in seq_along(idx)) {
    s   <- idx[ii]; rho <- post$rho[s]
    Kvv <- exp(-Dmat / rho); diag(Kvv) <- diag(Kvv) + 1e-8; Kvv_i <- solve(Kvv)
    L_Sigma <- diag(c(post$sigma_f[s, 1], post$sigma_f[s, 2])) %*% post$L_Omega[s, , ]
    for (g in seq_along(geom)) {
      G  <- geom[[g]]
      cd <- .cond_draw_icm(post$f_p[s, ], post$f_c[s, ], Kvv_i, G$Dgv, G$Dgg, rho, L_Sigma)
      eta_p <- post$b0_p[s] + cd$fp; eta_c <- post$b0_c[s] + cd$fc
      if (Pc > 0) {
        Xg    <- dat$X_grid[G$cells, , drop = FALSE]
        eta_p <- eta_p + as.numeric(Xg %*% post$beta_p[s, ])
        eta_c <- eta_c + as.numeric(Xg %*% post$beta_c[s, ])
      }
      pi_g <- plogis(eta_p)
      su   <- surv_cells("lognormal", eta_c, sigma_ln = post$sigma_ln[s], Tval = cfg$gpd_u)
      base[ii, g] <- weighted.mean(pi_g * su, G$pop)        # q_g with R = 1
    }
  }
  list(base = base, idx = idx)
}

log_msg("Per-LGA base draws (one pass)...")
pb  <- predict_lga_base(post, dat, geom)
base <- pb$base; idx <- pb$idx

# ---- per-draw GPD ratio at each threshold, then q_g(T) = base * R(T) ----
summ <- map_dfr(THR, function(T) {
  Rv <- vapply(idx, function(s) gpd_ratio(gpd_post$sigma[s], gpd_post$xi[s], CFG$gpd_u, T), 0)
  qg <- sweep(base, 1, Rv, "*")
  summarise_qg(qg, geom) |> mutate(threshold = T)
})
write_csv(summ, file.path(OUT, "q_g_multithreshold.csv"))
log_msg("Wrote q_g_multithreshold.csv")

# ---- check: T = 20,000 column reproduces q_g_summaries.csv ----
ref <- read_csv(file.path(OUT, "q_g_summaries.csv"), show_col_types = FALSE)
chk <- summ |> filter(threshold == 20000) |>
  inner_join(ref, by = c("NAME_1", "NAME_2"), suffix = c("", "_ref"))
cat(sprintf("CHECK vs q_g_summaries.csv (T=20k): max|d q_mean|=%.3e  max|d p_below|=%.3e\n",
            max(abs(chk$q_mean - chk$q_mean_ref)),
            max(abs(chk$p_below_0.1pct - chk$p_below_0.1pct_ref))))

# ---- headline numbers per threshold ----
hd <- summ |> group_by(threshold) |>
  summarise(max_per100k = max(q_mean) * 1e5,
            min_p_below = min(p_below_0.1pct),
            n_lga_p_below_lt1 = sum(p_below_0.1pct < 1),
            n_lga_cross_0.1pct = sum(q_mean >= CFG$ix_threshold), .groups = "drop")
cat("\nHEADLINE per threshold (all 78 LGAs):\n"); print(as.data.frame(hd))

# surveyed-only headline
surv_keys <- dat$d |> distinct(NAME_1, NAME_2)
hd_s <- summ |> inner_join(surv_keys, by = c("NAME_1","NAME_2")) |> group_by(threshold) |>
  summarise(n_surv = n(), max_per100k = max(q_mean)*1e5,
            min_p_below = min(p_below_0.1pct), .groups = "drop")
cat("\nHEADLINE per threshold (25 surveyed LGAs):\n"); print(as.data.frame(hd_s))
