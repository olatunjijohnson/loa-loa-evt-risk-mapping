# =============================================================================
# NGA_lga_fix_build.R — build a corrected survey file keyed on the SURVEY's own
# village identity (VillID), to test whether the LGA/coordinate labelling issue
# changes the headline results.
#
# Problem: villages are currently keyed on (NAME_1, NAME_2, Vill), i.e. on the
# GADM polygon each GPS point falls in. Several villages carry GPS points in two
# different LGAs -- four of them tens of km apart (data-entry errors) -- so those
# villages are split into two "villages" at two locations, giving 116 villages
# instead of the 110 the survey recorded (Emukah et al. 2018).
#
# Rule (deliberately minimal, to change only what is actually wrong):
#   For each surveyed village (VillID), keep the recorded points whose GADM LGA
#   matches the survey's own LGA label and discard the rest; if no point matches,
#   keep all of them (we cannot adjudicate). The village location is then the
#   mean of the kept rows -- the same averaging the pipeline already does, so
#   villages with ordinary within-village GPS scatter are left untouched.
#
# Writes data/NGA_lgafix.csv (same schema) -- originals untouched.
# =============================================================================
suppressMessages({ library(dplyr); library(readr); library(stringr) })

TS <- c("Abia","Anambra","Delta","Ebonyi","Imo")
d <- read_csv("data/NGA_cleaned_corrected.csv", show_col_types = FALSE) |>
  filter(NAME_1 %in% TS, !is.na(utmx), !is.na(utmy), !is.na(Loa_count))

norm <- function(x) str_squish(toupper(str_replace_all(x, "[^A-Za-z ]", " ")))

d <- d |> mutate(.match = norm(NAME_2) == norm(LGA))

# keep survey-consistent points where any exist; otherwise keep everything
keep <- d |> group_by(VillID) |> mutate(.any = any(.match)) |>
  filter(!.any | .match) |> ungroup()

vloc <- keep |> group_by(VillID) |>
  summarise(utmx_f = mean(utmx), utmy_f = mean(utmy),
            long_f = mean(long), lat_f = mean(lat),
            NAME_1_f = names(sort(table(NAME_1), decreasing = TRUE))[1],
            NAME_2_f = names(sort(table(NAME_2), decreasing = TRUE))[1],
            Vill_f   = names(sort(table(Vill),   decreasing = TRUE))[1],
            .groups = "drop")

old <- d |> group_by(NAME_1, NAME_2, Vill) |>
  summarise(ox = mean(utmx), oy = mean(utmy), n = n(), .groups = "drop")

fix <- d |> select(-.match) |> left_join(vloc, by = "VillID") |>
  mutate(utmx = utmx_f, utmy = utmy_f, long = long_f, lat = lat_f,
         NAME_1 = NAME_1_f, NAME_2 = NAME_2_f, Vill = Vill_f) |>
  select(-ends_with("_f"))

cat("rows dropped from location averaging:", nrow(d) - nrow(keep),
    "(their survey responses are retained; only the village coordinate changes)\n\n")

shift <- old |> left_join(vloc, by = character(0)) |> filter(FALSE)  # placeholder
cmp <- d |> distinct(VillID) |> left_join(vloc, by = "VillID") |>
  left_join(d |> group_by(VillID) |> summarise(ax = mean(utmx), ay = mean(utmy),
             Vill = first(Vill), LGA = first(LGA), .groups="drop"), by = "VillID") |>
  mutate(km_vs_allrow_mean = sqrt((utmx_f-ax)^2 + (utmy_f-ay)^2))
cat("village-location shift vs averaging ALL its rows (km):\n")
print(round(quantile(cmp$km_vs_allrow_mean, c(.5,.9,.95,1)), 3))
cat("\nvillages shifted more than 1 km by the correction:\n")
print(cmp |> filter(km_vs_allrow_mean > 1) |>
      transmute(VillID, Vill, surveyLGA = LGA, correctedLGA = NAME_2_f,
                shift_km = round(km_vs_allrow_mean,1)) |>
      arrange(desc(shift_km)) |> as.data.frame())

stopifnot(nrow(fix) == nrow(d))
cat("\nAFTER FIX: villages (N1,N2,Vill) =", nrow(distinct(fix, NAME_1, NAME_2, Vill)),
    "| VillID =", n_distinct(fix$VillID),
    "| GADM LGAs containing a village =", nrow(distinct(fix, NAME_1, NAME_2)), "\n")
cat("BEFORE   : villages =", nrow(old), "| GADM LGAs =", nrow(distinct(d, NAME_1, NAME_2)), "\n")

cv <- read_csv("data/covariates/covariates_villages_corrected.csv", show_col_types = FALSE)
miss <- distinct(fix, NAME_1, NAME_2, Vill) |> anti_join(cv, by = c("NAME_1","NAME_2","Vill"))
cat("village keys missing covariates:", nrow(miss), "\n")
if (nrow(miss)) print(as.data.frame(miss))

write_csv(fix, "data/NGA_lgafix.csv")
cat("\nWrote data/NGA_lgafix.csv\n")
