# =============================================================================
# NGA_report_figures.R — collaborator-facing maps, exploratory plots and a LaTeX
# macro file for the >20,000 mf/mL high-risk analysis. Uses LGA polygons from the
# ward shapefile (dissolved), a consistent blue(low)->red(high) palette, and per-cell
# surfaces computed with the SAME per-draw machinery as the per-LGA q_g.
#   Run AFTER NGA_run_highrisk.R.  Usage:  Rscript code/NGA_report_figures.R
# =============================================================================
source("code/NGA_common.R")
suppressMessages({ library(sf); library(ggplot2); library(dplyr); library(tidyr)
                   library(scales); library(patchwork) })
theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 12),
                plot.subtitle = element_text(size = 9, colour = "grey30"),
                legend.position = "right"))
# blue (low) -> red (high), per collaborator preference
BR <- c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B")
fill_br <- function(name, ..., na.value = "grey88") scale_fill_gradientn(colours = BR, name = name, na.value = na.value, ...)

OUT <- CFG$paths$out_dir; FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time()), "|", ..., "\n")

# -----------------------------------------------------------------------------
# 1. Data, fits, per-cell posterior-predictive surfaces (per-draw, like the pipeline)
# -----------------------------------------------------------------------------
log_msg("Loading fits + data...")
fr <- readRDS(file.path(OUT, "fits.rds")); dat <- fr$dat
post <- rstan::extract(fr$fits$lognormal)
gpd_post <- rstan::extract(fr$gpd)
d <- dat$d

log_msg("Per-cell surfaces (thinned per-draw)...")
Dmat <- dat$Dmat; Dgv <- dat$Dgv; Xg <- dat$X_grid
Ng <- nrow(dat$grid); Tval <- CFG$hr_threshold

# Per-cell posterior-mean surfaces from a model's posterior `pm`. tail in
# {weibull,lognormal,gamma,gpd}; "gpd" = lognormal body up to u then GPD beyond.
# Returns list(pi = prevalence, ix = proportion with mf > T), each length Ng.
cell_surface <- function(pm, tail, gpdp = NULL, thin = 10) {
  idx <- seq(1, length(pm$b0_p), by = thin); set.seed(1)
  pia <- numeric(Ng); ixa <- numeric(Ng)
  for (s in idx) {
    rho <- pm$rho[s]
    Kvv <- exp(-Dmat / rho); diag(Kvv) <- diag(Kvv) + 1e-8
    Kgv <- exp(-Dgv / rho); A <- Kgv %*% solve(Kvv)
    a   <- pmax(1 - rowSums(A * Kgv), 0)
    mp  <- as.numeric(A %*% pm$f_p[s, ]); mc <- as.numeric(A %*% pm$f_c[s, ])
    Ls  <- diag(c(pm$sigma_f[s, 1], pm$sigma_f[s, 2])) %*% pm$L_Omega[s, , ]
    res <- (matrix(rnorm(Ng * 2), Ng, 2) %*% t(Ls)) * sqrt(a)
    etap <- pm$b0_p[s] + mp + res[, 1] + as.numeric(Xg %*% pm$beta_p[s, ])
    etac <- pm$b0_c[s] + mc + res[, 2] + as.numeric(Xg %*% pm$beta_c[s, ])
    pis <- plogis(etap); pia <- pia + pis
    if (tail == "gpd") {
      R  <- gpd_ratio(gpdp$sigma[s], gpdp$xi[s], CFG$gpd_u, Tval)
      su <- plnorm(CFG$gpd_u, etac, pm$sigma_ln[s], lower.tail = FALSE)
      ixa <- ixa + pis * su * R
    } else {
      sv <- surv_cells(tail, etac, shape = pm$shape[s], sigma_ln = pm$sigma_ln[s], Tval = Tval)
      ixa <- ixa + pis * sv
    }
  }
  list(pi = pia / length(idx), ix = ixa / length(idx))
}

# headline surface = lognormal body + GPD tail; prevalence from the lognormal fit
cs <- cell_surface(post, "gpd", gpd_post)
grid <- dat$grid |> mutate(prevalence = cs$pi, hr_prop = cs$ix, hr_per100k = cs$ix * 1e5)

# -----------------------------------------------------------------------------
# 2. LGA / state / national polygons (dissolve wards; fix 4 name variants)
# -----------------------------------------------------------------------------
log_msg("Building boundaries...")
normk <- function(s, l) { k <- toupper(paste(s, l)); k <- gsub("-", " ", k)
  k <- gsub("[^A-Z ]", " ", k); gsub("\\s+", " ", trimws(k)) }
alias <- c("ABIA OBI NGWA" = "ABIA OBI NWGA", "IMO UNUIMO" = "IMO ONUIMO")  # grid -> shapefile
fix_key <- function(k) ifelse(k %in% names(alias), alias[k], k)

shp <- st_read("data/NGA_Ward_Boundaries_8721420816419582471/grid3_nga_boundary_vaccwards.shp",
               quiet = TRUE) |> st_make_valid()
states5 <- toupper(unique(dat$grid$NAME_1))
nation  <- shp |> group_by(statename) |> summarise(.groups = "drop")   # 37 states
state_poly <- nation[toupper(nation$statename) %in% states5, ]
shp5 <- shp[toupper(shp$statename) %in% states5, ]
lga_poly <- shp5 |> mutate(key = normk(statename, lganame)) |>
  group_by(key) |> summarise(.groups = "drop")

# -----------------------------------------------------------------------------
# 3. LGA-level tables: observed + predicted, joined to polygons
# -----------------------------------------------------------------------------
qs <- readr::read_csv(file.path(OUT, "q_g_summaries.csv"), show_col_types = FALSE) |>
  mutate(key = fix_key(normk(NAME_1, NAME_2)), hr_per100k = q_mean * 1e5)
obs_lga <- d |> group_by(NAME_1, NAME_2) |>
  summarise(obs_prev = mean(y_binary), n = n(), nvill = n_distinct(vid), .groups = "drop") |>
  mutate(key = fix_key(normk(NAME_1, NAME_2)))
pred_prev_lga <- grid |> group_by(NAME_1, NAME_2) |>
  summarise(prev = mean(prevalence), .groups = "drop") |>
  mutate(key = fix_key(normk(NAME_1, NAME_2)))

lga_obs  <- left_join(lga_poly, obs_lga, by = "key")
lga_prev <- left_join(lga_poly, pred_prev_lga, by = "key")
lga_hr   <- left_join(lga_poly, qs, by = "key")

# village lon/lat from individual records
vill <- d |> group_by(vid) |> summarise(long = mean(long), lat = mean(lat),
                                         prev = mean(y_binary), .groups = "drop")

dw <- median(diff(sort(unique(grid$long)))); th <- max(median(diff(sort(unique(grid$lat)))), .018)
tw <- max(dw, .018)
bb <- st_bbox(state_poly)
map_theme <- function() coord_sf(xlim = bb[c(1,3)], ylim = bb[c(2,4)], expand = TRUE)
lab_xy <- labs(x = NULL, y = NULL)

# =============================================================================
# FIGURE: study area
# =============================================================================
# Main panel zoomed to the five states (the study area is a small part of the
# country, so a whole-Nigeria view loses all detail), with a locator inset showing
# where that block sits nationally.
log_msg("Fig: study area (zoomed, with national locator inset)...")
f_main <- ggplot() +
  geom_sf(data = nation, fill = "grey96", colour = "grey80", linewidth = 0.15) +
  geom_sf(data = state_poly, aes(fill = statename), colour = "grey25", linewidth = 0.45) +
  geom_point(data = vill, aes(long, lat), size = 1.1, colour = "black", alpha = 0.75) +
  scale_fill_brewer(palette = "Set2", name = "Study state") +
  map_theme() + lab_xy +
  labs(title = "Study area: five states of south-eastern Nigeria",
       subtitle = "Survey villages (points) within Abia, Anambra, Delta, Ebonyi and Imo") +
  theme(legend.position = "right")
f_locator <- ggplot() +
  geom_sf(data = nation, fill = "grey93", colour = "grey70", linewidth = 0.12) +
  geom_sf(data = state_poly, fill = "#D6604D", colour = "grey25", linewidth = 0.15) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "grey55", linewidth = 0.4),
        plot.margin = margin(2, 2, 2, 2))
f_area <- f_main +
  inset_element(f_locator, left = 0.76, bottom = 0.02, right = 1.00, top = 0.30,
                align_to = "panel")
ggsave(file.path(FIG, "fig0_studyarea.png"), f_area, width = 8, height = 6.5, dpi = 200)

# =============================================================================
# FIGURE: exploratory data analysis
# =============================================================================
log_msg("Fig: EDA...")
cnt <- d$y_intensity[dat$pos_idx]
p_hist <- ggplot(data.frame(c = cnt), aes(c)) +
  geom_histogram(bins = 40, fill = "#4393C3", colour = "white") +
  geom_vline(xintercept = Tval, linetype = "dashed", colour = "#B2182B", linewidth = 1.1) +
  scale_x_log10(labels = comma) +
  labs(title = "A. Microfilarial density (positives)", x = "mf/mL (log scale)", y = "Individuals")
p_cdf <- ggplot(data.frame(c = cnt), aes(c)) + stat_ecdf(colour = "#2166AC", linewidth = 0.8) +
  geom_vline(xintercept = Tval, linetype = "dashed", colour = "#B2182B", linewidth = 1.1) +
  geom_vline(xintercept = max(cnt), linetype = "dotted", colour = "grey40", linewidth = 0.9) +
  annotate("text", x = max(cnt), y = 0.2, label = paste0("max ", comma(max(cnt))),
           hjust = -0.05, size = 3, colour = "grey30") +
  scale_x_log10(labels = comma) + scale_y_continuous(labels = percent) +
  labs(title = "B. Cumulative distribution", x = "mf/mL (log scale)", y = "Share of positives ≤ x")
p_prev <- ggplot(obs_lga, aes(reorder(NAME_2, obs_prev), obs_prev)) +
  geom_col(fill = "#4393C3") + coord_flip() + scale_y_continuous(labels = percent) +
  labs(title = "C. Observed prevalence by surveyed LGA", x = NULL, y = "Prevalence") +
  theme(axis.text.y = element_text(size = 6))
p_n <- ggplot(obs_lga, aes(n)) + geom_histogram(bins = 15, fill = "#92C5DE", colour = "white") +
  labs(title = "D. Individuals sampled per LGA", x = "n individuals", y = "LGAs")
f_eda <- (p_hist | p_cdf) / (p_prev | p_n) + plot_annotation(
  title = "Exploratory summary of the survey data",
  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(FIG, "fig_eda.png"), f_eda, width = 11, height = 9, dpi = 200)

# =============================================================================
# FIGURE 1: observed prevalence at LGA level (polygons; surveyed LGAs only)
# =============================================================================
log_msg("Fig 1: observed LGA prevalence...")
f1 <- ggplot(lga_obs) +
  geom_sf(aes(fill = obs_prev), colour = "grey70", linewidth = 0.15) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.5) +
  geom_point(data = vill, aes(long, lat), size = 0.5, colour = "black", alpha = 0.5) +
  fill_br("Observed\nprevalence", labels = percent_format(accuracy = 1)) +
  map_theme() + lab_xy +
  labs(title = "Observed Loa loa prevalence by LGA",
       subtitle = "Surveyed LGAs only (grey = not surveyed); points = survey villages")
ggsave(file.path(FIG, "fig1_obs_lga.png"), f1, width = 8, height = 6.5, dpi = 200)

# =============================================================================
# FIGURE 2: predicted prevalence — continuous grid + LGA aggregate (side by side)
# =============================================================================
log_msg("Fig 2: predicted prevalence (grid + LGA)...")
prev_lim <- range(c(grid$prevalence, lga_prev$prev), na.rm = TRUE)
g2a <- ggplot(grid) + geom_tile(aes(long, lat, fill = prevalence), width = tw, height = th) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.4, inherit.aes = FALSE) +
  fill_br("Prevalence", labels = percent_format(accuracy = 1), limits = prev_lim) +
  map_theme() + lab_xy + labs(title = "A. Continuous surface (2 km)")
g2b <- ggplot(lga_prev) + geom_sf(aes(fill = prev), colour = "grey70", linewidth = 0.15) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.4) +
  fill_br("Prevalence", labels = percent_format(accuracy = 1), limits = prev_lim, na.value = "white") +
  map_theme() + lab_xy + labs(title = "B. Aggregated to LGA")
f2 <- (g2a | g2b) + plot_layout(guides = "collect") +
  plot_annotation(title = "Predicted Loa loa prevalence",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(FIG, "fig2_prev.png"), f2, width = 12, height = 6, dpi = 200)

# =============================================================================
# FIGURES 3-4: predicted proportion with mf>20,000 — grid + LGA (side by side)
# =============================================================================
log_msg("Fig 3-4: high-risk proportion (grid + LGA)...")
hr_lim <- range(c(grid$hr_per100k, lga_hr$hr_per100k), na.rm = TRUE)
g3 <- ggplot(grid) + geom_tile(aes(long, lat, fill = hr_per100k), width = tw, height = th) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.4, inherit.aes = FALSE) +
  fill_br("per 100,000", trans = "log10", limits = hr_lim, labels = label_number(accuracy = .01)) +
  map_theme() + lab_xy + labs(title = "A. Continuous surface (2 km)")
g4 <- ggplot(lga_hr) + geom_sf(aes(fill = hr_per100k), colour = "grey70", linewidth = 0.15) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.4) +
  fill_br("per 100,000", trans = "log10", limits = hr_lim, labels = label_number(accuracy = .01), na.value = "white") +
  map_theme() + lab_xy + labs(title = "B. Aggregated to LGA")
f34 <- (g3 | g4) + plot_layout(guides = "collect") +
  plot_annotation(title = "Predicted proportion of people with mf > 20,000 mf/mL (per 100,000)",
                  subtitle = "Headline model: lognormal body + GPD tail. Note: 0.1% = 100 per 100,000 — all values are far below this.",
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.subtitle = element_text(size = 9, colour = "grey30")))
ggsave(file.path(FIG, "fig34_highrisk.png"), f34, width = 12, height = 6, dpi = 200)

# =============================================================================
# FIGURE: robustness across tail models (clear point-range vs threshold)
# =============================================================================
log_msg("Fig: robustness...")
ts <- readr::read_csv(file.path(OUT, "tail_sensitivity.csv"), show_col_types = FALSE) |>
  mutate(tail = recode(tail, weibull = "Weibull", lognormal = "Lognormal",
                       gamma = "Gamma", gpd = "Generalised Pareto (EVT)"))
ts_sum <- ts |> group_by(tail) |>
  summarise(lo = min(q_mean) * 1e5, md = median(q_mean) * 1e5,
            hi = max(q_mean) * 1e5, .groups = "drop")
thr <- CFG$ix_threshold * 1e5
f_rob <- ggplot(ts_sum, aes(y = reorder(tail, hi))) +
  geom_vline(xintercept = thr, linetype = "dashed", colour = "#B2182B") +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "#4393C3", linewidth = 1.2) +
  geom_point(aes(x = md), colour = "#2166AC", size = 2.5) +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000),
                labels = c("0.01","0.03","0.1","0.3","1","3","10","30","100","300","1,000"),
                limits = c(min(ts_sum$lo) * 0.5, thr * 1.5)) +
  labs(title = "Robustness of the high-risk burden across intensity-tail models",
       subtitle = "Across-LGA range (blue line) and median (dot) of the predicted >20,000 mf/mL burden per 100,000.\nThe red dashed line is the 0.1% threshold (100 per 100,000); even the most conservative tail stays far below it.",
       x = "Predicted people per 100,000 with mf > 20,000 (log scale)", y = NULL)
ggsave(file.path(FIG, "fig_robust.png"), f_rob, width = 9, height = 4.2, dpi = 200)

# =============================================================================
# FIGURE: ranked per-LGA caterpillar (predicted >20,000 burden, 95% CrI)
# =============================================================================
log_msg("Fig: ranked LGA caterpillar...")
qc <- qs |> arrange(q_mean) |> mutate(lab = paste0(NAME_2, " (", NAME_1, ")"),
        lo = q_lo * 1e5, mid = q_mean * 1e5, hi = q_hi * 1e5)
f_cat <- ggplot(qc, aes(x = mid, y = reorder(lab, mid))) +
  geom_vline(xintercept = thr, linetype = "dashed", colour = "#B2182B") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#92C5DE") +
  geom_point(colour = "#2166AC", size = 1.2) +
  scale_x_log10(labels = label_number(accuracy = .01)) +
  labs(title = "Predicted high-risk burden per LGA (ranked, 95% credible interval)",
       subtitle = "Red dashed line = 0.1% threshold (100 per 100,000); every LGA lies far to its left",
       x = "Predicted people per 100,000 with mf > 20,000", y = NULL) +
  theme(axis.text.y = element_text(size = 5))
ggsave(file.path(FIG, "fig_caterpillar.png"), f_cat, width = 8, height = 9, dpi = 200)

# LaTeX table: ALL surveyed LGAs (25), ranked by P(q_g<0.1%) ascending (lowest first).
top <- qs |> filter(key %in% obs_lga$key) |>
  arrange(p_below_0.1pct, desc(q_mean)) |>
  transmute(NAME_1, NAME_2, mean = sprintf("%.2f", q_mean * 1e5),
            cri = sprintf("%.2f--%.2f", q_lo * 1e5, q_hi * 1e5),
            pb = sprintf("%.4f", p_below_0.1pct))
tex <- c("\\begin{tabular}{llrrr}", "\\toprule",
  "State & LGA & Per 100{,}000 & 95\\% CrI & $P(q_g<0.1\\%)$ \\\\", "\\midrule",
  apply(top, 1, function(r) sprintf("%s & %s & %s & %s & %s \\\\",
        r["NAME_1"], r["NAME_2"], r["mean"], r["cri"], r["pb"])),
  "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(OUT, "top_lga_table.tex"))

# LaTeX table: LOO comparison of body distributions (basis for "lognormal best")
loo_f <- file.path(OUT, "validation", "F_loo_comparison.csv")
if (file.exists(loo_f)) {
  lc <- readr::read_csv(loo_f, show_col_types = FALSE) |> arrange(looic)
  ltex <- c("\\begin{tabular}{lrrr}", "\\toprule",
    "Body model & LOOIC & $p_{\\mathrm{loo}}$ & \\# bad $k$ \\\\", "\\midrule",
    apply(lc, 1, function(r) sprintf("%s & %.1f & %.1f & %s \\\\",
          tools::toTitleCase(r[["family"]]), as.numeric(r[["looic"]]), as.numeric(r[["p_loo"]]),
          if ("bad_k" %in% names(lc)) format(round(as.numeric(r[["bad_k"]]))) else "--")),
    "\\bottomrule", "\\end{tabular}")
  writeLines(ltex, file.path(OUT, "loo_table.tex"))
}

# remove stale figures from the previous report version
invisible(file.remove(Sys.glob(file.path(FIG, c("map[1-4]*.png", "chart[12]*.png")))))

# =============================================================================
# PARAMETER TABLE (param_table.tex)
# =============================================================================
log_msg("Parameter table...")
ci3 <- function(x) c(mean(x), quantile(x, c(.025, .975)))
prow <- function(lab, x, fmt = "%.2f") { v <- ci3(x)
  sprintf("%s & %s & %s--%s \\\\", lab, sprintf(fmt, v[1]), sprintf(fmt, v[2]), sprintf(fmt, v[3])) }
ptab <- c("\\begin{tabular}{lrr}", "\\toprule",
  "Parameter & Posterior mean & 95\\% CrI \\\\", "\\midrule",
  prow("Spatial range $\\rho$ (km)", post$rho, "%.1f"),
  prow("Prevalence field SD $\\sigma_{f,p}$", post$sigma_f[, 1]),
  prow("Intensity field SD $\\sigma_{f,c}$", post$sigma_f[, 2]),
  prow("Prevalence--intensity corr.\\ $\\rho_{pc}$", post$rho_pc),
  prow("Lognormal SD $\\sigma_{\\ln}$", post$sigma_ln),
  prow("GPD shape $\\xi$", gpd_post$xi),
  prow("GPD scale $\\sigma_u$ ($u=1{,}000$)", gpd_post$sigma, "%.0f"),
  "\\bottomrule", "\\end{tabular}")
writeLines(ptab, file.path(OUT, "param_table.tex"))

# =============================================================================
# FIGURE: horseshoe covariate-effect magnitudes (per-covariate spline-block norm)
# =============================================================================
log_msg("Fig: horseshoe covariate effects...")
meta <- dat$cov_meta
hs <- function(beta) do.call(rbind, lapply(unique(meta), function(cc) {
  cols <- which(meta == cc); nrm <- sqrt(rowSums(beta[, cols, drop = FALSE]^2))
  data.frame(cov = cc, m = mean(nrm), lo = quantile(nrm, .05), hi = quantile(nrm, .95)) }))
hsdf <- bind_rows(mutate(hs(post$beta_p), component = "Prevalence"),
                  mutate(hs(post$beta_c), component = "Intensity"))
f_hs <- ggplot(hsdf, aes(m, reorder(cov, m), colour = component)) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(.5)) +
  geom_point(position = position_dodge(.5)) +
  scale_colour_manual(values = c(Prevalence = "#2166AC", Intensity = "#B2182B"), name = NULL) +
  labs(title = "Covariate effect magnitudes under the horseshoe prior",
       subtitle = "Posterior magnitude of each covariate's spline-coefficient block (90% CrI); near zero = shrunk out",
       x = "Effect magnitude (standardised covariates)", y = NULL)
ggsave(file.path(FIG, "fig_horseshoe.png"), f_hs, width = 8, height = 5, dpi = 200)

# =============================================================================
# FIGURE: per-tail high-risk LGA maps (appendix), shared blue->red scale
# =============================================================================
log_msg("Fig: per-tail maps...")
TLEV <- c("Weibull", "Lognormal", "Gamma", "Generalised Pareto (EVT)")
ts_key <- ts |> mutate(key = fix_key(normk(NAME_1, NAME_2)), per100k = q_mean * 1e5,
                       tail = factor(tail, levels = TLEV))
# ONE shared log colour scale across panels (clean fixed breaks; tiny values floored).
scale_shared <- function() scale_fill_gradientn(colours = BR, name = "per 100,000",
  trans = "log10", limits = c(0.005, 7), breaks = c(0.01, 0.1, 1, 2, 5),
  labels = c("0.01", "0.1", "1", "2", "5"), oob = scales::squish, na.value = "grey88")
f_pertail <- lga_poly |> left_join(ts_key, by = "key") |> filter(!is.na(tail)) |>
  ggplot() +
  geom_sf(aes(fill = per100k), colour = "grey80", linewidth = 0.08) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.3) +
  facet_wrap(~tail) + scale_shared() + map_theme() + lab_xy +
  labs(title = "Predicted proportion exceeding 20,000 mf/mL, by LGA",
       subtitle = "Per 100,000 population, by intensity-tail model (shared log scale; blue = low, red = high)")
ggsave(file.path(FIG, "fig_pertail.png"), f_pertail, width = 11, height = 8.5, dpi = 200)

# =============================================================================
# FIGURE: per-tail high-risk CONTINUOUS (2 km grid) surfaces, per-panel log scales
# =============================================================================
log_msg("Fig: per-tail grid surfaces...")
postW2 <- rstan::extract(fr$fits$weibull); postG2 <- rstan::extract(fr$fits$gamma)
surf <- list(
  Weibull                   = cell_surface(postW2, "weibull")$ix,
  Lognormal                 = cell_surface(post,   "lognormal")$ix,
  Gamma                     = cell_surface(postG2, "gamma")$ix,
  `Generalised Pareto (EVT)` = cell_surface(post,  "gpd", gpd_post)$ix)
grid_long <- bind_rows(lapply(TLEV, function(tl)
  dat$grid |> transmute(long, lat, tail = tl, per100k = surf[[tl]] * 1e5))) |>
  mutate(tail = factor(tail, levels = TLEV))
f_pertail_grid <- ggplot(grid_long) +
  geom_tile(aes(long, lat, fill = per100k), width = tw, height = th) +
  geom_sf(data = state_poly, fill = NA, colour = "grey15", linewidth = 0.3, inherit.aes = FALSE) +
  facet_wrap(~tail) + scale_shared() + map_theme() + lab_xy +
  labs(title = "Predicted proportion exceeding 20,000 mf/mL, 2 km grid",
       subtitle = "Per 100,000 population, by intensity-tail model (shared log scale; blue = low, red = high)")
ggsave(file.path(FIG, "fig_pertail_grid.png"), f_pertail_grid, width = 11, height = 8.5, dpi = 200)

# =============================================================================
# FIGURE: tail diagnostics — MRL (u choice) + survival/extrapolation to 20,000
# =============================================================================
log_msg("Fig: tail diagnostics...")
us <- seq(200, 3500, 100)
mrl <- lapply(us, function(u) { e <- cnt[cnt > u] - u
  data.frame(u = u, me = mean(e), se = sd(e) / sqrt(length(e)), n = length(e)) }) |> bind_rows() |>
  filter(n >= 5)
p_mrl <- ggplot(mrl, aes(u, me)) +
  geom_ribbon(aes(ymin = me - 2 * se, ymax = me + 2 * se), alpha = .15, fill = "#2166AC") +
  geom_line(colour = "#2166AC") +
  geom_vline(xintercept = CFG$gpd_u, linetype = "dashed", colour = "#B2182B") +
  labs(title = "A. GPD threshold (mean residual life)",
       subtitle = "Roughly linear above a valid threshold; dashed = chosen u = 1,000",
       x = "threshold u (mf/mL)", y = "mean excess")
# representative model-implied conditional survival, plug-in posterior means
etac_v <- mean(post$b0_c) + colMeans(post$f_c) + as.numeric(dat$Xv %*% colMeans(post$beta_c))
sln <- mean(post$sigma_ln); u <- CFG$gpd_u; xi <- mean(gpd_post$xi); sg <- mean(gpd_post$sigma)
xs <- 10^seq(1, log10(25000), length = 250)
S_ln  <- sapply(xs, function(x) mean(plnorm(x, etac_v, sln, lower.tail = FALSE)))
Su    <- mean(plnorm(u, etac_v, sln, lower.tail = FALSE))
S_gpd <- ifelse(xs <= u, S_ln, Su * pmax(1 + xi * (xs - u) / sg, 0)^(-1 / xi))
# the other two body extrapolations, each from its OWN fit
postW <- rstan::extract(fr$fits$weibull); postG <- rstan::extract(fr$fits$gamma)
etacW <- mean(postW$b0_c) + colMeans(postW$f_c) + as.numeric(dat$Xv %*% colMeans(postW$beta_c))
etacG <- mean(postG$b0_c) + colMeans(postG$f_c) + as.numeric(dat$Xv %*% colMeans(postG$beta_c))
shW <- mean(postW$shape); shG <- mean(postG$shape)
S_w <- sapply(xs, function(x) mean(surv_cells("weibull", etacW, shape = shW, Tval = x)))
S_g <- sapply(xs, function(x) mean(surv_cells("gamma",   etacG, shape = shG, Tval = x)))
ecdf_pos <- ecdf(cnt)
sdf <- bind_rows(
  data.frame(x = xs, S = 1 - ecdf_pos(xs), model = "Empirical (positives)"),
  data.frame(x = xs, S = S_w,   model = "Weibull body"),
  data.frame(x = xs, S = S_g,   model = "Gamma body"),
  data.frame(x = xs, S = S_ln,  model = "Lognormal body"),
  data.frame(x = xs, S = S_gpd, model = "Lognormal + GPD tail")) |> filter(S > 0) |>
  mutate(model = factor(model, levels = c("Empirical (positives)", "Weibull body",
         "Gamma body", "Lognormal body", "Lognormal + GPD tail")))
p_surv <- ggplot(sdf, aes(x, S, colour = model)) + geom_line(linewidth = .8) +
  geom_vline(xintercept = max(cnt), linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 20000, linetype = "dashed", colour = "#B2182B") +
  scale_x_log10(labels = comma) + scale_y_log10(labels = scales::label_scientific()) +
  scale_colour_manual(values = c("Empirical (positives)" = "grey30",
        "Weibull body" = "#1B9E77", "Gamma body" = "#E6AB02",
        "Lognormal body" = "#2166AC", "Lognormal + GPD tail" = "#B2182B"), name = NULL) +
  labs(title = "B. Tail survival and extrapolation to 20,000 mf/mL",
       subtitle = "P(density > x | infected). Dotted = observed max (11,429); dashed = 20,000 threshold",
       x = "density x (mf/mL, log scale)", y = "P(density > x) (log)")
f_tail <- p_mrl + p_surv + plot_layout(widths = c(1, 1.2))
ggsave(file.path(FIG, "fig_tail.png"), f_tail, width = 12, height = 4.8, dpi = 200)

# =============================================================================
# LaTeX macro file with key numbers
# =============================================================================
log_msg("Writing LaTeX stats macros...")
lga_res <- readr::read_csv(file.path(OUT, "lga_results.csv"), show_col_types = FALSE)
val <- tryCatch(jsonlite::read_json(file.path(OUT, "validation/metrics.json")), error = function(e) NULL)
cv <- if (!is.null(val)) val$convergence[[1]] else NULL; pc <- if (!is.null(val)) val$prevalence_cv[[1]] else NULL
m <- c(
  sprintf("\\newcommand{\\Nind}{%s}", comma(nrow(d))),
  sprintf("\\newcommand{\\Npos}{%s}", comma(length(dat$pos_idx))),
  sprintf("\\newcommand{\\Nvill}{%d}", dat$V),
  sprintf("\\newcommand{\\Nlgasurv}{%d}", n_distinct(obs_lga$key)),
  sprintf("\\newcommand{\\Nlgapred}{%d}", nrow(lga_res)),
  sprintf("\\newcommand{\\Nstates}{%d}", length(states5)),
  sprintf("\\newcommand{\\prevoverall}{%.1f}", 100 * mean(d$y_binary)),
  sprintf("\\newcommand{\\maxcount}{%s}", comma(max(cnt))),
  sprintf("\\newcommand{\\medcount}{%s}", comma(round(median(cnt)))),
  sprintf("\\newcommand{\\pbelowlo}{%.3f}", min(lga_res$p_below)),     # headline = lognormal+GPD
  sprintf("\\newcommand{\\pbelowhi}{%.3f}", max(lga_res$p_below)),
  sprintf("\\newcommand{\\hrmax}{%.2f}", max(qs$hr_per100k)),
  sprintf("\\newcommand{\\maxrhat}{%.3f}", if (!is.null(cv)) cv$max_rhat else NA),
  sprintf("\\newcommand{\\miness}{%.0f}", if (!is.null(cv)) cv$min_bulk_ess else NA),
  sprintf("\\newcommand{\\ndiv}{%.0f}", if (!is.null(cv)) cv$n_divergent else NA),
  sprintf("\\newcommand{\\aucval}{%.2f}", if (!is.null(pc)) pc$auc else NA),
  sprintf("\\newcommand{\\moran}{%.2f}", if (!is.null(val)) val$morans_I else NA)
)
writeLines(m, file.path(OUT, "report_stats.tex"))
log_msg("Done. Figures + report_stats.tex in", OUT)
