# =============================================================================
# NGA_run_highrisk.R
# Answers the collaborator's question:
#   "What is the probability that <0.1% of unmeasured individuals in the survey
#    area have >20,000 mf/mL?"  -- reported BY LGA across the 5 states.
#
# Pipeline:
#   1. Load data (5 states), fit the JOINT coregionalized GP hurdle under
#      Weibull / lognormal / gamma, plus a pooled GPD (POT/EVT) tail.
#   2. Krige the latent fields to the 2 km grid.
#   3. Per cell, per draw: ix(x) = P(infected|x) * P(count > 20,000 | infected, x).
#   4. Per LGA, per draw: q_g = mean_cells ix(x); report P(q_g < 0.1%).
#   5. Tail-robustness band across the 4 distributions + LOO comparison.
#
# RUN ONLY AFTER NGA_validation.R confirms the model is adequate.
# Working directory must be the Nigeria/ folder.
# =============================================================================

source("code/NGA_common.R")
suppressMessages(library(loo))

# Optional overrides set by a driver (e.g. run_lite.R) BEFORE sourcing this file.
if (exists("RUN_OVERRIDE")) for (nm in names(RUN_OVERRIDE)) CFG[[nm]] <- RUN_OVERRIDE[[nm]]
if (exists("OUT_DIR_OVERRIDE")) CFG$paths$out_dir <- OUT_DIR_OVERRIDE

dir.create(CFG$paths$out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

# -----------------------------------------------------------------------------
# 1. Data + fits
# -----------------------------------------------------------------------------
log_msg("Loading data (5 states)...")
dat <- load_nga()
log_msg("N =", nrow(dat$d), "| villages =", dat$V,
        "| positives =", length(dat$pos_idx),
        "| max count =", max(dat$d$y_intensity))

log_msg("Compiling Stan models...")
models <- compile_models()

body_families <- c("weibull", "lognormal", "gamma")
fits <- list()
for (fam in body_families) {
  log_msg("Fitting hurdle —", fam, "...")
  fits[[fam]] <- fit_hurdle(dat, fam, models$hurdle)
}
log_msg("Fitting pooled GPD tail (u =", CFG$gpd_u, ")...")
fit_gpd_obj <- fit_gpd(dat, models$gpd)

saveRDS(list(fits = fits, gpd = fit_gpd_obj, dat = dat),
        file.path(CFG$paths$out_dir, "fits.rds"))

# -----------------------------------------------------------------------------
# 2-3. EXACT coregionalized per-LGA prediction (per-draw joint conditional draws)
# -----------------------------------------------------------------------------
log_msg("Precomputing per-LGA geometry (cross-distances)...")
geom <- precompute_lga_geometry(dat)

posts    <- lapply(fits, rstan::extract)
gpd_post <- rstan::extract(fit_gpd_obj)

log_msg("Exact per-LGA exceedance draws per tail...")
qg_list <- list()
for (tail in c(body_families, "gpd")) {
  log_msg("  tail:", tail)
  fam_post <- if (tail == "gpd") posts[["lognormal"]] else posts[[tail]]
  qg_list[[tail]] <- predict_lga_exact(
    fam_post, dat, geom, tail = tail, Tval = CFG$hr_threshold,
    gpd_post = if (tail == "gpd") gpd_post else NULL)
}

# -----------------------------------------------------------------------------
# 4. LGA-level answer + tail-robustness band
# -----------------------------------------------------------------------------
log_msg("Summarising to LGA...")
lga_by_tail <- imap_dfr(qg_list,
  function(qg, tail) summarise_qg(qg, geom) |> mutate(tail = tail))

# HEADLINE per LGA = lognormal body + GPD tail (the principled extrapolation under the
# best-fitting body); the across-tail min/max are kept as a robustness band.
lga_results <- lga_by_tail |>
  group_by(NAME_1, NAME_2) |>
  summarise(
    n_cells          = first(n_cells),
    p_below          = p_below_0.1pct[tail == "gpd"],   # headline: lognormal + GPD
    q_mean           = q_mean[tail == "gpd"],
    p_below_min      = min(p_below_0.1pct),             # robustness band across tails
    p_below_max      = max(p_below_0.1pct),
    .groups = "drop"
  ) |>
  arrange(p_below)

readr::write_csv(lga_by_tail, file.path(CFG$paths$out_dir, "tail_sensitivity.csv"))
readr::write_csv(lga_results, file.path(CFG$paths$out_dir, "lga_results.csv"))

# Per-LGA posterior summary of the HEADLINE q_g (lognormal body + GPD tail) for reporting
q_summaries <- summarise_qg(qg_list[["gpd"]], geom)
readr::write_csv(q_summaries, file.path(CFG$paths$out_dir, "q_g_summaries.csv"))

# -----------------------------------------------------------------------------
# 5. LOO comparison across body families (tail adequacy is judged in validation)
# -----------------------------------------------------------------------------
log_msg("LOO comparison (intensity component)...")
loo_list <- lapply(body_families, function(fam) {
  ll <- loo::extract_log_lik(fits[[fam]], "log_lik_int", merge_chains = FALSE)
  loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
})
names(loo_list) <- body_families
loo_cmp <- loo::loo_compare(loo_list)
write.csv(loo_cmp, file.path(CFG$paths$out_dir, "loo_comparison.csv"))

log_msg("Done. Results in", CFG$paths$out_dir)
print(lga_results, n = Inf)

# Render the report (md + tex + pdf)
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  log_msg("Rendering report (md + tex + pdf)...")
  # knit_root_dir = project root so the report's relative paths resolve; params$out
  # points the report at THIS run's results dir (full run or subset).
  for (fmt in c("md_document", "latex_document", "pdf_document")) {
    try(rmarkdown::render("code/report.Rmd", output_format = fmt,
                          output_dir = CFG$paths$out_dir,
                          knit_root_dir = getwd(),
                          params = list(out = CFG$paths$out_dir), quiet = TRUE), silent = TRUE)
  }
}
