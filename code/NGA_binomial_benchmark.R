# =============================================================================
# NGA_binomial_benchmark.R — model-free, one-sided binomial upper limits on the
# >20,000 mf/mL high-risk proportion, per surveyed LGA. A transparent benchmark to
# sit alongside the model (Peter's request). For x exceedances in n people the
# one-sided 95% upper confidence limit is the Clopper-Pearson upper bound (the
# "rule of three", ~3/n, when x = 0). Contrasted with the model's 95% upper limit.
#   Run AFTER NGA_run_highrisk.R + NGA_multithreshold.R.
#   Usage:  Rscript code/NGA_binomial_benchmark.R
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(dplyr); library(readr) })
OUT <- CFG$paths$out_dir; T <- 20000L
ub <- function(x, n, conf = 0.95) if (x == 0) 1 - (1 - conf)^(1/n) else qbeta(conf, x + 1, n - x)

dat <- load_nga(); d <- dat$d
N <- nrow(d)

# pooled (whole study area) one-sided 95% upper limit at each threshold
for (Tt in c(10000, 15000, 20000)) {
  x <- sum(d$y_intensity > Tt)
  cat(sprintf("POOLED >%d: %d/%d exceed -> 95%% upper = %.1f per 100k\n", Tt, x, N, ub(x, N) * 1e5))
}

# per-LGA: people tested (n), exceedances (x), naive upper limit, model 95% upper (q_hi)
lga <- d |> group_by(NAME_1, NAME_2) |>
  summarise(n = n(), x = sum(y_intensity > T), .groups = "drop") |>
  mutate(naive_up = vapply(seq_len(n()), function(i) ub(x[i], n[i]), 0) * 1e5)

mt <- read_csv(file.path(OUT, "q_g_multithreshold.csv"), show_col_types = FALSE) |>
  filter(threshold == T) |> transmute(NAME_1, NAME_2, model_up = q_hi * 1e5)
lga <- lga |> inner_join(mt, by = c("NAME_1", "NAME_2")) |> arrange(desc(naive_up))

cat(sprintf("\nPer-LGA (%d surveyed): n median %d (range %d-%d)\n",
            nrow(lga), median(lga$n), min(lga$n), max(lga$n)))
cat(sprintf("Naive 95%% upper per 100k: median %.0f (range %.0f-%.0f); #below 100/100k bar: %d/%d\n",
            median(lga$naive_up), min(lga$naive_up), max(lga$naive_up),
            sum(lga$naive_up < 100), nrow(lga)))
cat(sprintf("Model 95%% upper per 100k: median %.1f (max %.1f); #below 100/100k bar: %d/%d\n",
            median(lga$model_up), max(lga$model_up), sum(lga$model_up < 100), nrow(lga)))

# LaTeX table
fmt <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
rows <- apply(lga, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s \\\\",
  r["NAME_1"], r["NAME_2"], fmt(as.numeric(r["n"])), as.character(as.integer(r["x"])),
  fmt(as.numeric(r["naive_up"])), sprintf("%.1f", as.numeric(r["model_up"]))))
tex <- c("\\begin{tabular}{llrrrr}", "\\toprule",
  "State & LGA & People & Above & Simple 95\\% & Model 95\\% \\\\",
  " &  & tested & 20{,}000 & upper limit & upper limit \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(OUT, "binomial_table.tex"))
cat("\nWrote binomial_table.tex\n")
