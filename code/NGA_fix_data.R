# =============================================================================
# NGA_fix_data.R — correct two mis-geocoded villages and write corrected inputs.
#   * "Edo"  rows are really Delta / Ethiope East  (village Okpara Inland)
#   * "Benue" rows are really Ebonyi / Ohaukwu     (village Ndiakparata Ugwuagba)
# The survey LGA column is the trusted label; coordinates were verified to fall
# inside the correct LGA polygons, so they are kept. The two villages were absent
# from the GEE covariate extract, so they are given their nearest grid cell's
# covariates. Originals are left untouched; corrected copies are written.
# =============================================================================
suppressMessages({ library(dplyr); library(readr) })
source("code/NGA_common.R")   # for COVARS

raw <- read_csv("data/Emmanuel_cleaned_NGA_data.csv", show_col_types = FALSE)

# --- relabel from the trusted LGA column ---------------------------------------
fix <- raw |>
  mutate(
    NAME_2 = case_when(NAME_1 %in% c("Edo", "Benue") ~ LGA, TRUE ~ NAME_2),
    NAME_1 = case_when(NAME_1 == "Edo" ~ "Delta",
                       NAME_1 == "Benue" ~ "Ebonyi",
                       TRUE ~ NAME_1),
    NAME_0 = "Nigeria")
stopifnot(all(c("Delta","Ebonyi") %in% fix$NAME_1), !any(fix$NAME_1 %in% c("Edo","Benue")))
write_csv(fix, "data/NGA_cleaned_corrected.csv")

# --- new villages that need covariates (now in target states, absent from cov) -
ts <- c("Abia","Anambra","Delta","Ebonyi","Imo")
fix5 <- fix |> filter(NAME_1 %in% ts, !is.na(utmx), !is.na(utmy), !is.na(Loa_count))
cvill <- read_csv("data/covariates/covariates_villages.csv", show_col_types = FALSE)
grid  <- read_csv("data/covariates/covariates_grid.csv", show_col_types = FALSE)

newv <- fix5 |>
  group_by(NAME_0, NAME_1, NAME_2, Vill) |>
  summarise(long = mean(long), lat = mean(lat), utmx = mean(utmx), utmy = mean(utmy),
            .groups = "drop") |>
  anti_join(cvill, by = c("NAME_1","NAME_2","Vill"))
cat("New villages needing covariates:\n"); print(newv[, c("NAME_1","NAME_2","Vill")])

# nearest grid cell (Euclidean in km) -> copy its covariates
gx <- as.matrix(grid[, c("utmx","utmy")])
nn_cov <- function(ux, uy) {
  i <- which.min((gx[,1]-ux)^2 + (gx[,2]-uy)^2)
  grid[i, COVARS]
}
newrows <- newv |> rowwise() |>
  mutate(nn_cov(utmx, utmy)) |> ungroup()

# append to covariate-village file (fill non-covariate columns as needed)
add <- bind_rows(cvill[0, ], newrows)          # align columns to cvill schema
for (cc in COVARS) add[[cc]] <- newrows[[cc]]
add$NAME_0 <- newrows$NAME_0; add$NAME_1 <- newrows$NAME_1
add$NAME_2 <- newrows$NAME_2; add$Vill <- newrows$Vill
add$long <- newrows$long; add$lat <- newrows$lat
add$utmx <- newrows$utmx; add$utmy <- newrows$utmy
cvill_fixed <- bind_rows(cvill, add)
stopifnot(!any(is.na(cvill_fixed[COVARS])))
write_csv(cvill_fixed, "data/covariates/covariates_villages_corrected.csv")

cat("\nWrote data/NGA_cleaned_corrected.csv (", nrow(fix), "rows ) and",
    "covariates_villages_corrected.csv (", nrow(cvill_fixed), "villages )\n")
cat("Target-state individuals:", nrow(fix5),
    "| villages:", n_distinct(fix5$Vill, fix5$NAME_2, fix5$NAME_1),
    "| surveyed LGAs:", n_distinct(paste(fix5$NAME_1, fix5$NAME_2)),
    "| max count:", max(fix5$Loa_count), "\n")
