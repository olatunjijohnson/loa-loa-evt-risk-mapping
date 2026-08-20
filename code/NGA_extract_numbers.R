# =============================================================================
# NGA_extract_numbers.R — pull every number the manuscript quotes out of a
# results directory, so the text and tables can be updated from one place.
#   LOA_OUT=results/highrisk_20000_fixed Rscript code/NGA_extract_numbers.R
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(dplyr); library(readr) })
OUT <- CFG$paths$out_dir
cat("=== results dir:", OUT, "| data:", CFG$paths$data, "===\n\n")

fr <- readRDS(file.path(OUT, "fits.rds")); dat <- fr$dat; fit <- fr$fits$lognormal
d <- dat$d

cat("--- STUDY DESCRIPTION ---\n")
cat("individuals            :", nrow(d), "\n")
cat("villages               :", dat$V, "\n")
cat("surveyed LGAs          :", nrow(distinct(d, NAME_1, NAME_2)), "\n")
cat("positives (count>0)    :", length(dat$pos_idx),
    sprintf(" (%.2f%%)\n", 100*length(dat$pos_idx)/nrow(d)))
ci <- d$y_intensity[dat$pos_idx]
cat("median / max density   :", median(ci), "/", max(ci), "\n")
cat("exceedances over u     :", sum(ci > CFG$gpd_u), "\n")
bl <- d |> group_by(NAME_1, NAME_2) |> summarise(n=n(), prev=mean(y_binary), .groups="drop")
cat("LGAs with observed prev > 20% :", sum(bl$prev > .2), "->",
    paste(bl$NAME_2[bl$prev > .2], collapse=", "), "\n")
bv <- d |> count(vid); bg <- d |> count(NAME_1, NAME_2)
cat("per community: range", min(bv$n), "-", max(bv$n), "median", median(bv$n),
    "IQR", paste(quantile(bv$n, c(.25,.75)), collapse="-"), "\n")
cat("per LGA      : range", min(bg$n), "-", max(bg$n), "median", median(bg$n),
    "IQR", paste(quantile(bg$n, c(.25,.75)), collapse="-"), "\n\n")

cat("--- TABLE 1: PARAMETERS (mean, 95% CrI) ---\n")
sm <- rstan::summary(fit)$summary
gs <- rstan::summary(fr$gpd)$summary
row <- function(lbl, m, digits=2) {
  if (!m %in% rownames(sm) && !m %in% rownames(gs)) { cat(sprintf("%-28s NOT FOUND\n", lbl)); return(invisible()) }
  r <- if (m %in% rownames(sm)) sm[m,] else gs[m,]
  cat(sprintf("%-28s %s  (%s -- %s)\n", lbl, round(r["mean"],digits),
              round(r["2.5%"],digits), round(r["97.5%"],digits)))
}
row("Spatial range rho (km)", "rho", 1)
row("Prevalence field SD",    "sigma_f[1]")
row("Intensity field SD",     "sigma_f[2]")
om <- grep("^L_Omega\\[2,1\\]$|^Omega\\[1,2\\]$|^rho_pc$|^omega$", rownames(sm), value=TRUE)
if (length(om)) { for (o in om) row(paste("cross-corr", o), o) } else {
  po <- rstan::extract(fit)
  if (!is.null(po$L_Omega)) {
    L <- po$L_Omega; Om <- L[,2,1]*L[,1,1] + L[,2,2]*0   # Omega_21 = L21*L11
    Om <- apply(L, 1, function(M) (M %*% t(M))[1,2])
    cat(sprintf("%-28s %s  (%s -- %s)\n", "Prev-intensity corr",
        round(mean(Om),2), round(quantile(Om,.025),2), round(quantile(Om,.975),2)))
  }
}
row("Lognormal SD sigma_ln",  "sigma_ln")
row("GPD shape xi",           "xi")
row("GPD scale sigma_u",      "sigma", 0)

cat("\n--- CONVERGENCE ---\n")
cat("max Rhat      :", round(max(sm[,"Rhat"], na.rm=TRUE), 4), "\n")
cat("min n_eff     :", round(min(sm[,"n_eff"], na.rm=TRUE)), "\n")
nd <- rstan::get_num_divergent(fit); ni <- sum(sapply(rstan::get_sampler_params(fit, inc_warmup=FALSE), nrow))
cat("divergences   :", nd, sprintf("of %d post-warmup (%.1f%%)\n", ni, 100*nd/ni))
cat("min E-BFMI    :", round(min(rstan::get_bfmi(fit)), 3), "\n")
cat("max treedepth hits:", rstan::get_num_max_treedepth(fit), "\n")

cat("\n--- HEADLINE q_g ---\n")
qs <- read_csv(file.path(OUT, "q_g_summaries.csv"), show_col_types=FALSE)
cat("prediction LGAs        :", nrow(qs), "\n")
cat("max burden per 100,000 :", round(max(qs$q_mean)*1e5, 2),
    sprintf("(%.4f%%)\n", 100*max(qs$q_mean)))
cat("P(q<0.1%) range        :", round(min(qs$p_below_0.1pct),4), "-",
    round(max(qs$p_below_0.1pct),4), "\n")
cat("all LGAs > 0.99        :", all(qs$p_below_0.1pct > 0.99), "\n")
cat("highest-burden LGA     :", with(qs[which.max(qs$q_mean),], paste(NAME_1, NAME_2)), "\n")

cat("\n--- TABLE 2: SURVEYED LGAs (ranked by P then burden) ---\n")
surv <- distinct(d, NAME_1, NAME_2)
t2 <- qs |> semi_join(surv, by=c("NAME_1","NAME_2")) |>
  arrange(p_below_0.1pct, desc(q_mean)) |>
  transmute(NAME_1, NAME_2,
            per100k = sprintf("%.2f", q_mean*1e5),
            CrI = sprintf("%.2f--%.2f", q_lo*1e5, q_hi*1e5),
            P = sprintf("%.4f", p_below_0.1pct))
cat("rows:", nrow(t2), "\n")
for (i in seq_len(nrow(t2))) cat(sprintf("%s & %s & %s & %s & %s \\\\\n",
    t2$NAME_1[i], t2$NAME_2[i], t2$per100k[i], t2$CrI[i], t2$P[i]))
cat("\nworst surveyed P:", min(as.numeric(t2$P)), "| max surveyed burden:",
    max(as.numeric(t2$per100k)), "\n")

mt <- file.path(OUT, "q_g_multithreshold.csv")
if (file.exists(mt)) {
  cat("\n--- MULTITHRESHOLD ---\n")
  m <- read_csv(mt, show_col_types=FALSE)
  nm <- intersect(c("p_below_0.1pct","p_below"), names(m))[1]
  for (Tv in sort(unique(m$T))) {
    s <- m[m$T == Tv, ]
    ss <- s |> semi_join(surv, by=c("NAME_1","NAME_2"))
    cat(sprintf("T=%6d  max per100k=%7.2f  min P(all)=%.4f  min P(surveyed)=%.4f\n",
        Tv, max(s$q_mean)*1e5, min(s[[nm]]), min(ss[[nm]])))
  }
}
ts <- file.path(OUT, "tail_sensitivity.csv")
if (file.exists(ts)) {
  cat("\n--- TAIL SENSITIVITY (across-LGA range per 100,000) ---\n")
  print(read_csv(ts, show_col_types=FALSE) |> group_by(tail) |>
        summarise(lo=min(q_mean)*1e5, med=median(q_mean)*1e5, hi=max(q_mean)*1e5,
                  minP=min(p_below_0.1pct), .groups="drop") |> as.data.frame(), digits=3)
}
