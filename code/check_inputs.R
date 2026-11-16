# Read-only preflight for the corrected manuscript workflow; does not fit models.
covariates <- c("elevation", "slope", "lst_day", "ndvi", "rain_total",
                "rain_seasonality", "treecover", "water_occurrence", "population")
inputs <- list()
inputs[[Sys.getenv("LOA_DATA", "data/NGA_lgafix.csv")]] <-
  c("NAME_1", "NAME_2", "Vill", "utmx", "utmy", "Loa_count")
inputs[["data/covariates/covariates_villages_corrected.csv"]] <-
  c("NAME_1", "NAME_2", "Vill", covariates)
inputs[["data/covariates/covariates_grid.csv"]] <-
  c("NAME_1", "NAME_2", "utmx", "utmy", "long", "lat", covariates)
ok <- TRUE
for (path in names(inputs)) {
  if (!file.exists(path)) {
    cat("MISSING input:", path, "\n"); ok <- FALSE; next
  }
  fields <- names(read.csv(path, nrows = 0, check.names = FALSE))
  missing <- setdiff(inputs[[path]], fields)
  if (length(missing)) {
    cat("MISSING columns in", path, ":", paste(missing, collapse = ", "), "\n")
    ok <- FALSE
  } else cat("OK input header:", path, "\n")
}
shape <- "data/NGA_Ward_Boundaries_8721420816419582471/grid3_nga_boundary_vaccwards"
for (extension in c("shp", "shx", "dbf", "prj")) {
  path <- paste0(shape, ".", extension)
  if (!file.exists(path)) {
    cat("MISSING map input:", path, "\n"); ok <- FALSE
  }
}
packages <- c("tidyverse", "rstan", "loo", "fields", "pROC", "scoringRules",
              "sf", "patchwork", "jsonlite", "rmarkdown")
installed <- rownames(installed.packages())
missing <- setdiff(packages, installed)
if (length(missing)) {
  cat("MISSING R packages:", paste(missing, collapse = ", "), "\n"); ok <- FALSE
}
for (package in intersect(packages, installed)) {
  problem <- tryCatch({ loadNamespace(package); NULL }, error = conditionMessage)
  if (!is.null(problem)) {
    cat("FAILED TO LOAD", package, ":", problem, "\n"); ok <- FALSE
  }
}
if ("patchwork" %in% installed && "ggplot2" %in% installed &&
    packageVersion("ggplot2") >= "4.0.0" && packageVersion("patchwork") < "1.3.0") {
  cat("INCOMPATIBLE: ggplot2 >= 4 needs patchwork >= 1.3.0.\n"); ok <- FALSE
}
cat("Checked file presence, headers and package loading only; no data or models were changed.\n")
quit(status = if (ok) 0 else 1)
