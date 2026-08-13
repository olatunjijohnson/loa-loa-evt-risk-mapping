# =============================================================================
# NGA_manuscript_figures_rev.R — revised manuscript figures (co-author review).
#   (1) EDA figure: two panels only — the intensity ECDF (A) and observed
#       prevalence by surveyed LGA (B). The earlier histogram and the
#       "individuals per LGA" histogram were dropped; the sample-size ranges
#       that replace the latter are printed to the console for use in the text.
#   (2) Tail-robustness figure: a main-text version restricted to the lognormal
#       body and the headline lognormal+GPD tail, plus a supplementary version
#       retaining all four body/tail models (Weibull, gamma, lognormal,
#       lognormal+GPD).
# Reads only the cleaned survey CSV and tail_sensitivity.csv, so it re-runs in
# seconds without reloading the fitted Stan objects.
#
# Usage:  Rscript code/NGA_manuscript_figures_rev.R
#   Set LOA_MS_FIGDIR=/path/to/manuscript/figures to mirror the figures into a
#   manuscript directory as well as results/highrisk_20000/figures.
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr)
                   library(scales); library(patchwork) })
theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 12),
                plot.subtitle = element_text(size = 9, colour = "grey30"),
                legend.position = "right"))

OUT <- CFG$paths$out_dir; FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MSFIG <- Sys.getenv("LOA_MS_FIGDIR", unset = "")   # optional manuscript figures dir
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

# Save to results/, and mirror into the manuscript directory when one is set.
save_fig <- function(name, plot, ...) {
  ggsave(file.path(FIG, name), plot, ...)
  if (nzchar(MSFIG) && dir.exists(MSFIG)) ggsave(file.path(MSFIG, name), plot, ...)
}

# -----------------------------------------------------------------------------
# Individual records (same filtering / village keying as load_nga())
# -----------------------------------------------------------------------------
d <- readr::read_csv(CFG$paths$data, show_col_types = FALSE) |>
  filter(NAME_1 %in% CFG$target_states) |>
  filter(!is.na(utmx), !is.na(utmy), !is.na(Loa_count)) |>
  mutate(y_binary = as.integer(Loa_count > 0),
         y_intensity = as.numeric(Loa_count)) |>
  unite("vill_key", all_of(CFG$village_key), remove = FALSE, sep = " | ")

Tval <- CFG$hr_threshold
cnt  <- d$y_intensity[d$y_binary == 1]

by_vill <- d |> count(vill_key, name = "n")
by_lga  <- d |> count(NAME_1, NAME_2, name = "n")

log_msg("individuals:", nrow(d), "| villages:", nrow(by_vill), "| surveyed LGAs:", nrow(by_lga))
log_msg("per community: range", min(by_vill$n), "-", max(by_vill$n),
        "| median", median(by_vill$n), "| IQR", paste(quantile(by_vill$n, c(.25, .75)), collapse = "-"))
log_msg("per LGA:       range", min(by_lga$n), "-", max(by_lga$n),
        "| median", median(by_lga$n), "| IQR", paste(quantile(by_lga$n, c(.25, .75)), collapse = "-"))

# -----------------------------------------------------------------------------
# (1) Fig 2 — EDA with two panels: ECDF (A) and observed LGA prevalence (B)
# -----------------------------------------------------------------------------
log_msg("Fig: EDA (two panels)...")
p_cdf <- ggplot(data.frame(c = cnt), aes(c)) +
  stat_ecdf(colour = "#2166AC", linewidth = 0.8) +
  geom_vline(xintercept = Tval, linetype = "dashed", colour = "#B2182B", linewidth = 1.1) +
  geom_vline(xintercept = max(cnt), linetype = "dotted", colour = "grey40", linewidth = 0.9) +
  annotate("text", x = max(cnt), y = 0.2, label = paste0("max ", comma(max(cnt))),
           hjust = 1.05, size = 3, colour = "grey30") +
  annotate("text", x = Tval, y = 0.45, label = "20,000 mf/mL", angle = 90,
           vjust = -0.4, size = 3, colour = "#B2182B") +
  scale_x_log10(labels = comma, expand = expansion(mult = c(0.02, 0.10))) +
  scale_y_continuous(labels = percent) +
  labs(title = "A. Cumulative distribution", x = "mf/mL (log scale)", y = "Share of positives ≤ x")

obs_lga <- d |> group_by(NAME_1, NAME_2) |>
  summarise(obs_prev = mean(y_binary), n = n(), .groups = "drop")
p_prev <- ggplot(obs_lga, aes(reorder(NAME_2, obs_prev), obs_prev)) +
  geom_col(fill = "#4393C3") + coord_flip() + scale_y_continuous(labels = percent) +
  labs(title = "B. Observed prevalence by surveyed LGA", x = NULL, y = "Prevalence") +
  theme(axis.text.y = element_text(size = 7))

f_eda <- (p_cdf | p_prev) + plot_annotation(
  title = "Exploratory summary of the survey data",
  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig("fig_eda.png", f_eda, width = 11, height = 4.8, dpi = 200)

# -----------------------------------------------------------------------------
# (2) Fig 7 — robustness across tail models: two-model (main) and four (S1)
# -----------------------------------------------------------------------------
log_msg("Fig: robustness...")
lab <- c(weibull = "Weibull", lognormal = "Lognormal",
         gamma = "Gamma", gpd = "Lognormal + GPD tail")
ts <- readr::read_csv(file.path(OUT, "tail_sensitivity.csv"), show_col_types = FALSE) |>
  mutate(tail = recode(tail, !!!lab))
ts_sum <- ts |> group_by(tail) |>
  summarise(lo = min(q_mean) * 1e5, md = median(q_mean) * 1e5,
            hi = max(q_mean) * 1e5, .groups = "drop")
thr <- CFG$ix_threshold * 1e5
print(as.data.frame(ts_sum))

rob_plot <- function(dd, subtitle) {
  ggplot(dd, aes(y = reorder(tail, hi))) +
    geom_vline(xintercept = thr, linetype = "dashed", colour = "#B2182B") +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "#4393C3", linewidth = 1.2) +
    geom_point(aes(x = md), colour = "#2166AC", size = 2.5) +
    scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000),
                  labels = c("0.01","0.03","0.1","0.3","1","3","10","30","100","300","1,000"),
                  limits = c(min(dd$lo) * 0.5, thr * 1.5)) +
    labs(title = "Robustness of the high-risk burden across intensity-tail models",
         subtitle = subtitle, y = NULL,
         x = "Predicted people per 100,000 with mf > 20,000 (log scale)")
}
sub_main <- paste0("Across-LGA range (blue line) and median (dot) of the predicted >20,000 mf/mL burden per 100,000.\n",
                   "The red dashed line is the 0.1% threshold (100 per 100,000); the heavier GPD tail stays far below it.")
sub_all  <- paste0("Across-LGA range (blue line) and median (dot) of the predicted >20,000 mf/mL burden per 100,000.\n",
                   "The red dashed line is the 0.1% threshold (100 per 100,000); all four tail models stay far below it.")

# main text: lognormal body vs lognormal + GPD tail
f_rob <- rob_plot(filter(ts_sum, tail %in% c("Lognormal", "Lognormal + GPD tail")), sub_main)
save_fig("fig_robust.png", f_rob, width = 9, height = 3.2, dpi = 200)

# supplement: all four body/tail models
f_rob_all <- rob_plot(ts_sum, sub_all)
save_fig("fig_robust_all.png", f_rob_all, width = 9, height = 4.2, dpi = 200)

log_msg("Done.")
