# =============================================================================
# NGA_lga_fix_run.R — does correcting the village/LGA labelling change the
# headline result?
#
# Runs the headline model (lognormal body + pooled GPD tail) TWICE, through
# identical code: once on the ORIGINAL data (116 villages, GADM-keyed) and once
# on the CORRECTED data (110 villages, survey-keyed; NGA_lga_fix_build.R), then
# compares the per-LGA q_g and P(q_g < 0.1%).
#
# One patch is required. build_covariates() places natural-spline interior knots
# at quantiles of the village-level covariate values. water_occurrence is zero at
# most villages, so at 116 villages its interior knots are (0, 0, 0.25) and at
# 110 villages they collapse to (0, 0, 0), which splines::ns() rejects. We
# therefore place interior knots at quantiles of the DISTINCT values and drop any
# that still coincide with a boundary, falling back to a linear term if none
# survive. Both runs use this same patched basis, so the comparison is like for
# like.  Writes results/lgafix_compare/.
# =============================================================================
source("code/NGA_common.R")
OUT <- "results/lgafix_compare"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

# ---- patched covariate basis (robust interior knots) ------------------------
build_covariates <- function(villages, grid, cfg = CFG) {
  covs  <- if (!is.null(cfg$cov_subset)) cfg$cov_subset else COVARS
  cvill <- readr::read_csv(cfg$paths$cov_vill, show_col_types = FALSE)
  prep  <- function(D) dplyr::mutate(D,
            population = log1p(ifelse(is.na(population), 0, population)))
  cvill <- prep(cvill)
  vjoin <- villages |>
    left_join(cvill |> select(NAME_1, NAME_2, Vill, all_of(covs)),
              by = c("NAME_1", "NAME_2", "Vill"))
  if (any(is.na(vjoin[covs]))) stop("Missing covariate values for some villages after join.")
  gridc <- prep(grid)
  mu <- vapply(covs, function(c) mean(vjoin[[c]]), numeric(1))
  sg <- vapply(covs, function(c) sd(vjoin[[c]]),   numeric(1))

  Bv <- list(); Bg <- list(); meta <- character()
  for (c in covs) {
    xv <- (vjoin[[c]] - mu[c]) / sg[c]
    xg <- (gridc[[c]] - mu[c]) / sg[c]
    bk <- range(xv)
    k  <- unique(as.numeric(quantile(unique(xv), seq(0, 1, length.out = cfg$spline_df + 1)[-c(1, cfg$spline_df + 1)])))
    k  <- k[k > bk[1] & k < bk[2]]
    if (length(k)) {
      b <- splines::ns(xv, knots = k, Boundary.knots = bk)
    } else {
      b <- matrix(xv, ncol = 1); attr(b, "linear") <- TRUE
      message("  linear fallback for covariate: ", c)
    }
    if (!is.null(attr(b, "linear"))) {
      Bv[[c]] <- b; Bg[[c]] <- matrix(xg, ncol = 1)
    } else {
      Bv[[c]] <- b; Bg[[c]] <- predict(b, xg)
    }
    if (length(k) < cfg$spline_df - 1)
      message("  reduced basis for ", c, ": ", ncol(Bv[[c]]), " columns (interior knots kept: ", length(k), ")")
    meta <- c(meta, rep(c, ncol(Bv[[c]])))
  }
  Xv <- do.call(cbind, Bv); Xg <- do.call(cbind, Bg)
  colnames(Xv) <- colnames(Xg) <- paste0(meta, ".", ave(meta, meta, FUN = seq_along))
  list(Xv = Xv, Xg = Xg, meta = meta, mu = mu, sd = sg)
}

run_one <- function(tag, datafile) {
  log_msg("=============== ", tag, " (", datafile, ") ===============")
  CFG$paths$data <<- datafile
  dat <- load_nga()
  log_msg(tag, ": N =", nrow(dat$d), "| villages =", dat$V,
          "| positives =", length(dat$pos_idx), "| covariate cols =", ncol(dat$X))
  models <- compile_models()
  log_msg(tag, ": fitting lognormal hurdle...")
  fit_ln <- fit_hurdle(dat, "lognormal", models$hurdle)
  log_msg(tag, ": fitting GPD tail...")
  fit_gp <- fit_gpd(dat, models$gpd)

  sm <- rstan::summary(fit_ln)$summary
  gs <- rstan::summary(fit_gp)$summary
  log_msg(tag, ": max Rhat =", round(max(sm[,"Rhat"], na.rm=TRUE), 4),
          "| min n_eff =", round(min(sm[,"n_eff"], na.rm=TRUE)))
  pr <- rbind(sm[intersect(c("rho","sigma_f[1]","sigma_f[2]","sigma_ln"), rownames(sm)),
                 c("mean","2.5%","97.5%"), drop=FALSE],
              gs[intersect(c("xi","sigma"), rownames(gs)), c("mean","2.5%","97.5%"), drop=FALSE])
  cat("\n", tag, "key parameters:\n"); print(round(pr, 3))

  geom <- precompute_lga_geometry(dat)
  post <- rstan::extract(fit_ln); gpp <- rstan::extract(fit_gp)
  log_msg(tag, ": exact per-LGA prediction...")
  qg <- predict_lga_exact(post, dat, geom, tail = "gpd", Tval = CFG$hr_threshold, gpd_post = gpp)
  qs <- summarise_qg(qg, geom)
  readr::write_csv(qs, file.path(OUT, paste0("q_g_", tag, ".csv")))
  saveRDS(list(par = pr, V = dat$V, npos = length(dat$pos_idx), ncov = ncol(dat$X)),
          file.path(OUT, paste0("meta_", tag, ".rds")))
  log_msg(tag, ": max burden/100k =", round(max(qs$q_mean)*1e5, 3),
          "| P(q<0.1%) range", round(min(qs$p_below_0.1pct),4), "-", round(max(qs$p_below_0.1pct),4))
  qs
}

a <- run_one("original",  "data/NGA_cleaned_corrected.csv")
b <- run_one("corrected", "data/NGA_lgafix.csv")

cmp <- full_join(a |> select(NAME_1, NAME_2, q_o = q_mean, p_o = p_below_0.1pct),
                 b |> select(NAME_1, NAME_2, q_c = q_mean, p_c = p_below_0.1pct),
                 by = c("NAME_1","NAME_2")) |>
  mutate(per100k_o = q_o*1e5, per100k_c = q_c*1e5,
         ratio = per100k_c/per100k_o, dp = p_c - p_o) |>
  arrange(desc(per100k_c))
readr::write_csv(cmp, file.path(OUT, "comparison.csv"))

cat("\n\n================ COMPARISON ================\n")
cat("LGAs: original", nrow(a), "| corrected", nrow(b), "\n")
cat("max burden per 100,000 : original", round(max(a$q_mean)*1e5,3),
    "| corrected", round(max(b$q_mean)*1e5,3), "\n")
cat("min P(q_g<0.1%)        : original", round(min(a$p_below_0.1pct),4),
    "| corrected", round(min(b$p_below_0.1pct),4), "\n")
cat("all LGAs > 0.99        : original", all(a$p_below_0.1pct>0.99),
    "| corrected", all(b$p_below_0.1pct>0.99), "\n")
cat("\nTop 15 LGAs after correction:\n")
print(cmp |> head(15) |> transmute(NAME_1, NAME_2,
      before = round(per100k_o,2), after = round(per100k_c,2),
      ratio = round(ratio,2), P_before = round(p_o,4), P_after = round(p_c,4)) |> as.data.frame())
cat("\nburden ratio (corrected/original) quantiles:\n")
print(round(quantile(cmp$ratio, c(0,.25,.5,.75,1), na.rm=TRUE), 3))
log_msg("Done.")
