# Loa loa safety classification in Nigeria

Analysis code accompanying **A joint geostatistical extreme-value model to enable
safe ivermectin mass drug administration in Loa loa-endemic areas**, by Olatunji
Johnson, Emmanuel Emukah, Lindsay J. Rakers, Jenna Coalson, Charles Nzelu,
Bertram E. B. Nwoke, Katie Gass and Peter J. Diggle. The manuscript is being
prepared for submission; no publication DOI is available yet.

The study retrospectively analyses the **2016** LoaScope survey in Abia, Anambra,
Delta, Ebonyi and Imo. The manuscript analysis uses 10,605 people in 110 villages,
with 654 positive measurements. The corrected village locations fall in 24
administrative polygons; predictions cover 78 LGAs. These are historical,
model-dependent estimates, not a current operational risk map.

The headline analysis combines a bivariate geostatistical model of infection
prevalence and positive intensity with a separately fitted, pooled generalised
Pareto tail above 1,000 mf/mL. It estimates `q_g`, the population-weighted fraction
exceeding 20,000 mf/mL in LGA `g`, and `P(q_g < 0.001)`. The manuscript compares
this probability with the proposed provisional 99% criterion. A joint spatial-tail
extension is implemented separately; it is not the headline fit.

## Data access

**The survey data, covariate files, fitted models and generated results are not
included.** A fresh clone therefore cannot reproduce the analysis by itself.
The individual-level records contain potentially identifying location information.
Their access arrangements remain those described in the manuscript: requests are
handled by the study custodians under the original consent and ethics terms.
Public release is being discussed with the survey team; it has not been authorised
or performed here.

The original survey is [Emukah et al. (2018)](https://doi.org/10.4269/ajtmh.18-0163).
Related earlier modelling is [Johnson et al. (2022)](https://doi.org/10.1371/journal.pntd.0010189).

## Required inputs

Run commands from the repository root. The corrected manuscript workflow expects:

| File | Purpose and required fields |
|---|---|
| `data/NGA_lgafix.csv` | Corrected survey records: `NAME_1`, `NAME_2`, `Vill`, `utmx`, `utmy`, `Loa_count` |
| `data/covariates/covariates_villages_corrected.csv` | Village keys `NAME_1`, `NAME_2`, `Vill`, plus the nine covariates below |
| `data/covariates/covariates_grid.csv` | Prediction grid with `NAME_1`, `NAME_2`, `utmx`, `utmy`, `long`, `lat`, plus the nine covariates |
| `data/NGA_Ward_Boundaries_8721420816419582471/grid3_nga_boundary_vaccwards.shp` | Map boundaries; include matching `.dbf`, `.shx` and `.prj` files. Attributes include `statename` and `lganame` |

Coordinates `utmx` and `utmy` are in **kilometres**, not metres. The covariates are
`elevation`, `slope`, `lst_day`, `ndvi`, `rain_total`, `rain_seasonality`,
`treecover`, `water_occurrence` and `population`. Prediction aggregation uses
population weights. The map scripts dissolve **GRID3 Nigeria ward polygons** into
LGA and state boundaries; GADM labels in the survey serve a different purpose.
The exact boundary release and reuse terms should accompany any redistributed
boundary files; none are bundled here.

`code/NGA_gee_extract.py` documents environmental extraction with Google Earth
Engine. It requires the study locations, an authenticated account, your own Cloud
project ID, and Python packages `earthengine-api` and `pandas`. It does not recreate
restricted survey locations. `data/prediction_NGA_states.csv` is the fallback
prediction grid for runs without covariates; those runs are not the manuscript analysis.

## Software

Use R >= 4.1 and a working C++ toolchain for `rstan`. Required R packages are:

```r
install.packages(c("tidyverse", "rstan", "loo", "fields", "pROC",
                   "scoringRules", "sf", "patchwork", "jsonlite", "rmarkdown"))
```

The spatial package `sf` needs its platform's GDAL/GEOS/PROJ libraries.
Use `patchwork` >= 1.3.0 with ggplot2 4.x. The optional
`run_prevalence_cv_metrics.R` also uses `PRROC`. HTML reports use Pandoc;
PDF reports need a working LaTeX installation. This repository does not yet
contain a complete frozen software environment, so record `sessionInfo()` for
each reproduction and allow for Monte Carlo and package-version differences.

## Reproduce the corrected headline analysis

First obtain authorised copies of the inputs and check their presence:

```bash
Rscript code/check_inputs.R
```

If starting from the earlier cleaned survey, run `Rscript code/NGA_lga_fix_build.R`
once to create `data/NGA_lgafix.csv`. This needs the earlier
`data/NGA_cleaned_corrected.csv`, including `VillID`, the survey `LGA`, `long`,
`lat` and the fields listed above. It retains participants while correcting the
coordinate used to represent each village. Review the correction rule before
applying it to a different survey.

Then run:

```bash
# Spatial-block and random-village cross-validation; inspect diagnostics.
bash code/run_lgafix_validation.sh

# Main fits, threshold comparisons, maps and tables.
bash code/run_lgafix_full.sh

# Final main/supplementary EDA and robustness figure variants.
LOA_DATA=data/NGA_lgafix.csv LOA_OUT=results/highrisk_20000_fixed \
  Rscript code/NGA_manuscript_figures_rev.R
```

Validation performs repeated MCMC fits and can take several hours. Successful
execution is not an automatic declaration of model adequacy: inspect convergence,
calibration and tail diagnostics. The main fit uses four chains, 2,000 iterations
per chain and 1,000 warm-up iterations; cross-validation uses shorter fits.
These settings and seed 1234 are in `code/NGA_common.R`.

**Avoid mixing the two result directories.** Corrected manuscript headline results
belong in `results/highrisk_20000_fixed/`. The older `results/highrisk_20000/`
directory and `run_corrected.sh` predate the village identity correction. Individual
R scripts retain the earlier default path for compatibility, so specify both
overrides when running them directly:

```bash
export LOA_DATA=data/NGA_lgafix.csv
export LOA_OUT=results/highrisk_20000_fixed
Rscript code/NGA_shape_stability.R
Rscript code/NGA_prior_sensitivity_xi.R
Rscript code/NGA_village_cv.R
Rscript code/NGA_binomial_benchmark.R
```

## Find the results

Paths below are relative to `$LOA_OUT`:

| Output | Contents |
|---|---|
| `q_g_summaries.csv` | Headline per-LGA posterior burdens, intervals and classification probabilities |
| `q_g_multithreshold.csv` | Results at 10,000, 15,000 and 20,000 mf/mL |
| `tail_sensitivity.csv` | Burdens under alternative positive-density/tail models |
| `validation/metrics.json` | Spatial-block CV and convergence summaries |
| `validation/village_cv_metrics.json` | Random-village CV summaries |
| `validation/intensity_cv_metrics.csv` | Held-out intensity diagnostics, saved by the updated validation script for future plot regeneration |
| `validation/C_loo_pit.pdf`, `.eps`, `.png` | Vector and 600-dpi exports of the spatial-block CV PIT plot; legacy filename retained |
| `figures/` | Maps and summary plots |
| `fits.rds` | Fitted Stan models and supporting data; keep access-controlled |

Figures can be mirrored into a manuscript directory by setting `LOA_MS_FIGDIR`
for scripts that support it. `NGA_extract_numbers.R` provides a numerical audit;
it does not substitute for inspecting the diagnostics.

## Spatial-tail extension and interpretation

`NGA_joint_tail.R` fits pooled and spatial-scale variants jointly with the hurdle
model. `NGA_joint_tail_compare.R` compares their per-LGA summaries with the
headline fit in the same `$LOA_OUT`. `NGA_spatial_tail_prototype.R` is an earlier
two-stage prototype. These are distinct analyses: do not substitute their output
for the headline model or describe the separately fitted headline as the joint
extension. The extension's likelihood description and prior documentation remain
under manuscript review.

The corrected 20,000 mf/mL analysis gives a minimum classification probability of
0.998 across 78 LGAs. At 10,000 mf/mL the minimum is 0.983, so the same 99%
criterion is not met everywhere. The prior-sensitivity script examines conditional
tail survival; it does not by itself establish robustness of every LGA's full
classification probability. These distinctions matter when interpreting outputs.

## License

The code is available under the [MIT license](LICENSE). This license does not
grant access to, or permission to redistribute, third-party survey or boundary data.
