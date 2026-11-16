# =============================================================================
# NGA_joint_tail_compare.R   (Plan 3 / Task 5: new vs previous results)
# -----------------------------------------------------------------------------
# Compares the spatially-linked-tail per-LGA answer against the committed
# headline (lognormal body + pooled separate-GPD tail) in q_g_summaries.csv.
# Run AFTER NGA_joint_tail.R (full) has written joint_qg_summaries.csv.
# Output: results/highrisk_20000/spatial_tail/joint_vs_previous.{csv,md}
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(tidyr) })
ROOT <- Sys.getenv("LOA_ROOT", unset = getwd())  # run from the repository root, or set LOA_ROOT
setwd(ROOT)
OUTD <- Sys.getenv("LOA_OUT", "results/highrisk_20000")
OUT  <- file.path(OUTD, "spatial_tail")

prev  <- read_csv(file.path(OUTD, "q_g_summaries.csv"), show_col_types = FALSE) |>
  transmute(NAME_1, NAME_2,
            q_prev = q_mean, p_prev = p_below_0.1pct)        # committed headline (separate GPD)
joint <- read_csv(file.path(OUT, "joint_qg_summaries.csv"), show_col_types = FALSE)

sp <- joint |> filter(variant == "spatial_tail") |>
  transmute(NAME_1, NAME_2, q_sp = q_mean, p_sp = p_below_0.1pct)
po <- joint |> filter(variant == "pooled_tail_joint") |>
  transmute(NAME_1, NAME_2, q_po = q_mean, p_po = p_below_0.1pct)

cmp <- prev |> left_join(sp, by = c("NAME_1","NAME_2")) |>
  left_join(po, by = c("NAME_1","NAME_2")) |>
  mutate(
    burden_prev_100k = q_prev * 1e5,
    burden_spatial_100k = q_sp * 1e5,
    d_burden_100k = (q_sp - q_prev) * 1e5,
    d_p_below = p_sp - p_prev
  ) |>
  arrange(desc(abs(d_burden_100k)))

write_csv(cmp, file.path(OUT, "joint_vs_previous.csv"))

con <- file(file.path(OUT, "joint_vs_previous.md"), "w")
wl <- function(...) cat(..., "\n", file = con, sep = "")
wl("# Spatially-linked tail vs committed headline\n")
wl("Comparison of the per-LGA high-risk answer under three tail treatments:\n")
wl("- **prev**: committed headline (lognormal body + pooled GPD fitted *separately*).")
wl("- **pooled_joint**: same pooled tail, but fitted *jointly* with the hurdle (sanity check).")
wl("- **spatial**: Peter's tail, sigma_u(x)=exp(alpha+beta*f_c(x)), fitted jointly.\n")
wl(sprintf("LGAs compared: %d\n", nrow(cmp)))
wl("## Headline answer P(q_g < 0.1%)\n")
wl(sprintf("- prev    : min %.4f  median %.4f  max %.4f", min(cmp$p_prev), median(cmp$p_prev), max(cmp$p_prev)))
wl(sprintf("- pooled_joint : min %.4f  median %.4f  max %.4f", min(cmp$p_po), median(cmp$p_po), max(cmp$p_po)))
wl(sprintf("- spatial : min %.4f  median %.4f  max %.4f\n", min(cmp$p_sp), median(cmp$p_sp), max(cmp$p_sp)))
wl("## Worst-LGA burden (per 100,000)\n")
wl(sprintf("- prev    : %.3f", max(cmp$burden_prev_100k)))
wl(sprintf("- spatial : %.3f\n", max(cmp$burden_spatial_100k, na.rm = TRUE)))
wl(sprintf("Largest absolute burden change (spatial - prev): %.3f per 100k\n", max(abs(cmp$d_burden_100k), na.rm = TRUE)))
wl("## 10 LGAs that move most (spatial vs prev)\n")
top <- head(cmp, 10) |>
  transmute(LGA = paste(NAME_1, NAME_2, sep = " / "),
            burden_prev = round(burden_prev_100k, 3),
            burden_spatial = round(burden_spatial_100k, 3),
            d_burden = round(d_burden_100k, 3),
            p_prev = round(p_prev, 4), p_spatial = round(p_sp, 4))
wl("| LGA | burden_prev | burden_spatial | Δburden | P_prev | P_spatial |")
wl("|---|---|---|---|---|---|")
for (i in seq_len(nrow(top))) wl(sprintf("| %s | %.3f | %.3f | %+.3f | %.4f | %.4f |",
  top$LGA[i], top$burden_prev[i], top$burden_spatial[i], top$d_burden[i], top$p_prev[i], top$p_spatial[i]))
close(con)
cat(readLines(file.path(OUT, "joint_vs_previous.md")), sep = "\n")
cat("\nWrote", file.path(OUT, "joint_vs_previous.csv"), "and .md\n")
