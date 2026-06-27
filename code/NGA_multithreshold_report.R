# =============================================================================
# NGA_multithreshold_report.R — table + faceted map for the 10k/15k/20k mf/mL
# exceedance analysis. Reads q_g_multithreshold.csv (from NGA_multithreshold.R).
#   Usage:  Rscript code/NGA_multithreshold_report.R
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(sf); library(ggplot2); library(dplyr); library(tidyr)
                   library(scales); library(readr); library(stringr) })
OUT <- CFG$paths$out_dir; FIG <- file.path(OUT, "figures")
theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 12),
                legend.position = "right"))
BR <- c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B")

mt  <- read_csv(file.path(OUT, "q_g_multithreshold.csv"), show_col_types = FALSE)
dat <- load_nga()

# ---- surveyed LGAs (the 25), ranked by P(q_g<0.1%) at the strictest test (10k) ----
normk <- function(s, l) { k <- toupper(paste(s, l)); k <- gsub("-", " ", k)
  k <- gsub("[^A-Z ]", " ", k); gsub("\\s+", " ", trimws(k)) }
alias <- c("ABIA OBI NGWA" = "ABIA OBI NWGA", "IMO UNUIMO" = "IMO ONUIMO")
fix_key <- function(k) ifelse(k %in% names(alias), alias[k], k)
surv_keys <- dat$d |> distinct(NAME_1, NAME_2)

w <- mt |> inner_join(surv_keys, by = c("NAME_1","NAME_2")) |>
  mutate(per100k = q_mean * 1e5) |>
  select(NAME_1, NAME_2, threshold, per100k, p_below_0.1pct) |>
  pivot_wider(names_from = threshold, values_from = c(per100k, p_below_0.1pct))
ord <- order(w$`p_below_0.1pct_10000`, -w$`per100k_10000`)
w <- w[ord, ]

fnum <- function(x) sprintf("%.2f", x)
fp   <- function(x) sprintf("%.4f", x)
rows <- apply(w, 1, function(r) sprintf(
  "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
  r["NAME_1"], r["NAME_2"],
  fnum(as.numeric(r["per100k_10000"])), fnum(as.numeric(r["per100k_15000"])), fnum(as.numeric(r["per100k_20000"])),
  fp(as.numeric(r["p_below_0.1pct_10000"])), fp(as.numeric(r["p_below_0.1pct_15000"])), fp(as.numeric(r["p_below_0.1pct_20000"]))))
tex <- c(
  "\\begin{tabular}{ll rrr rrr}", "\\toprule",
  " &  & \\multicolumn{3}{c}{Burden per 100{,}000} & \\multicolumn{3}{c}{$P(q_g<0.1\\%)$} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "State & LGA & 10k & 15k & 20k & 10k & 15k & 20k \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(OUT, "multithreshold_table.tex"))
cat("Wrote multithreshold_table.tex (", nrow(w), "surveyed LGAs)\n")

# ---- faceted LGA map: 10k / 15k / 20k, shared log scale, white = no prediction ----
shp <- st_read("data/NGA_Ward_Boundaries_8721420816419582471/grid3_nga_boundary_vaccwards.shp",
               quiet = TRUE) |> st_make_valid()
states5 <- toupper(unique(dat$grid$NAME_1))
nation  <- shp |> group_by(statename) |> summarise(.groups = "drop")
state_poly <- nation[toupper(nation$statename) %in% states5, ]
shp5 <- shp[toupper(shp$statename) %in% states5, ]
lga_poly <- shp5 |> mutate(key = normk(statename, lganame)) |> group_by(key) |> summarise(.groups = "drop")

flab <- function(T) paste0("mf > ", format(T, big.mark = ","), " mf/mL")
flev <- flab(c(10000, 15000, 20000))
mt2  <- mt |> mutate(key = fix_key(normk(NAME_1, NAME_2)), per100k = q_mean * 1e5)
lim  <- range(mt2$per100k)
# expand every polygon across all 3 thresholds so unpredicted LGAs render white in
# each panel (rather than collapsing into a stray NA facet).
mp <- bind_rows(lapply(c(10000, 15000, 20000), function(T) {
  lga_poly |> left_join(filter(mt2, threshold == T), by = "key") |> mutate(facet = flab(T))
})) |> mutate(facet = factor(facet, levels = flev))
bb  <- st_bbox(state_poly)
f_mt <- ggplot(mp) +
  geom_sf(aes(fill = per100k), colour = "grey70", linewidth = 0.12) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.35) +
  facet_wrap(~ facet, nrow = 1) +
  scale_fill_gradientn(colours = BR, name = "per 100,000", trans = "log10",
                       limits = lim, labels = label_number(accuracy = 0.01), na.value = "white") +
  coord_sf(xlim = bb[c(1,3)], ylim = bb[c(2,4)], expand = TRUE) +
  labs(title = "Predicted high-risk burden at three mf/mL thresholds",
       subtitle = "LGA population-weighted mean of P(mf > threshold), per 100,000. Shared log scale; white = no prediction.") +
  theme(axis.text = element_text(size = 6), strip.text = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_multithreshold.png"), f_mt, width = 12, height = 4.6, dpi = 200)
cat("Wrote fig_multithreshold.png\n")
