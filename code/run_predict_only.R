# =============================================================================
# run_predict_only.R — recompute the per-LGA answer from the EXISTING fits, with
# the two updates: population-weighted LGA means and a lognormal+GPD headline.
# Reuses results/highrisk_20000/fits.rds (model unchanged) so no refit is needed.
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(dplyr); library(purrr); library(readr) })
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")
OUT <- CFG$paths$out_dir

fr <- readRDS(file.path(OUT, "fits.rds")); dat <- fr$dat
log_msg("Per-LGA geometry (population-weighted)...")
geom <- precompute_lga_geometry(dat)
posts <- lapply(fr$fits, rstan::extract); gpd_post <- rstan::extract(fr$gpd)

body_families <- c("weibull", "lognormal", "gamma")
qg_list <- list()
for (tail in c(body_families, "gpd")) {
  log_msg("  tail:", tail)
  fam_post <- if (tail == "gpd") posts[["lognormal"]] else posts[[tail]]
  qg_list[[tail]] <- predict_lga_exact(fam_post, dat, geom, tail = tail,
    Tval = CFG$hr_threshold, gpd_post = if (tail == "gpd") gpd_post else NULL)
}

lga_by_tail <- imap_dfr(qg_list, function(qg, tail) summarise_qg(qg, geom) |> mutate(tail = tail))
lga_results <- lga_by_tail |>
  group_by(NAME_1, NAME_2) |>
  summarise(n_cells = first(n_cells),
            p_below = p_below_0.1pct[tail == "gpd"], q_mean = q_mean[tail == "gpd"],
            p_below_min = min(p_below_0.1pct), p_below_max = max(p_below_0.1pct),
            .groups = "drop") |> arrange(p_below)
write_csv(lga_by_tail, file.path(OUT, "tail_sensitivity.csv"))
write_csv(lga_results, file.path(OUT, "lga_results.csv"))
write_csv(summarise_qg(qg_list[["gpd"]], geom), file.path(OUT, "q_g_summaries.csv"))
log_msg("Done. Headline = lognormal+GPD, population-weighted. P(<0.1%) range:",
        paste(round(range(lga_results$p_below), 4), collapse = "-"))
