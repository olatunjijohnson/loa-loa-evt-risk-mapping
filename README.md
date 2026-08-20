# Loa loa high-intensity risk mapping with extreme-value theory

Code for a Bayesian **joint geostatistical and extreme-value** model that predicts
where dangerously high *Loa loa* microfilarial loads (>20,000 mf/mL) are appreciably
common, to support the safety of ivermectin mass drug administration (MDA) for
onchocerciasis and lymphatic filariasis in loiasis co-endemic areas. The analysis is
applied to a survey of south-eastern Nigeria (Abia, Anambra, Delta, Ebonyi, Imo).

The model
1. maps infection **prevalence** (logistic geostatistical regression),
2. maps microfilarial **intensity** among positives (lognormal body), the two fields
   linked by a **bivariate Gaussian process** (intrinsic coregionalisation), and
3. characterises the **upper tail** of intensity with a peaks-over-threshold
   **generalised Pareto distribution (GPD)** — necessary because the threshold of
   concern (20,000 mf/mL) lies beyond the largest observed density.

For each Local Government Area (LGA) it reports `P(q_g < 0.1%)`, the posterior
probability that fewer than 0.1% of people exceed 20,000 mf/mL, with robustness
checks across alternative tail models, lower thresholds, and a model-free benchmark.
It also implements a **spatially-structured tail** extension,
`σ_u(x) = exp(α + β·f_c(x))`, that lets tail heaviness vary in space.

## Paper

> [Authors]. *Predicting where few people carry dangerously high Loa loa
> microfilarial loads: a Bayesian joint geostatistical and extreme-value analysis to
> support safe ivermectin mass treatment in south-eastern Nigeria.* (in preparation /
> under review).
> **DOI: _to be added on publication._**

Closely related prior work using geostatistics for loiasis MDA safety assessment:
Johnson O, Giorgi E, Fronterrè C, et al. *Geostatistical modelling enables efficient
safety assessment for mass drug administration with ivermectin in Loa loa endemic
areas...* PLOS Negl Trop Dis. 2022;16(2):e0010189.
<https://doi.org/10.1371/journal.pntd.0010189>

## Repository layout

```
code/
  NGA_common.R               Model definitions: Stan hurdle + ICM GP, GPD tail,
                             kriging, exact per-LGA prediction, helpers (CFG, load_nga)
  NGA_fix_data.R             Data-relabelling step (two misgeocoded villages)
  NGA_gee_extract.py         Google Earth Engine extraction of environmental covariates
  NGA_validation.R           Adequacy gate: convergence, calibration, spatial block CV,
                             tail (mean-residual-life) diagnostics
  NGA_run_highrisk.R         Main run: fit, krige, per-LGA q_g, tail-robustness, LOO
  NGA_report_figures.R       All maps and figures
  NGA_manuscript_figures_rev.R  Revised manuscript figures: two-panel EDA figure and
                             the tail-robustness figure in main-text (lognormal vs
                             lognormal+GPD) and supplementary (all four models) forms
  NGA_multithreshold.R       Burden / P(q_g<0.1%) at 10k, 15k, 20k mf/mL
  NGA_multithreshold_report.R  Multi-threshold table + faceted map
  NGA_binomial_benchmark.R   Model-free per-LGA "rule of three" upper limits
  NGA_shape_stability.R      GPD threshold-stability scan over u = 500..1500 mf/mL
  NGA_prior_sensitivity_xi.R Sensitivity of the tail shape xi to its prior
  NGA_village_cv.R           Random village-fold CV (local interpolation skill),
                             the counterpart to the harsh spatial block CV
  NGA_lga_fix_build.R        Rebuild the survey file keyed on the SURVEY's own village
                             identity (VillID), correcting villages whose records carry
                             GPS points in two different LGAs
  NGA_lga_fix_run.R          Headline model fitted twice through identical code, on the
                             original and the corrected input, with a per-LGA comparison
  NGA_extract_numbers.R      Pull every number the manuscript quotes from a results dir
  run_lgafix_full.sh         Driver: full pipeline on the corrected input
  run_lgafix_validation.sh   Driver: validation + village CV on the corrected input
  NGA_spatial_tail_prototype.R  Two-stage spatial-GPD prototype (σ_u(x)=exp(α+β f_c))
  NGA_joint_tail.R           Joint hurdle+spatial-GPD fit (β-free vs β=0) + per-LGA q_g
  NGA_joint_tail_compare.R   Spatial-tail vs pooled-tail per-LGA comparison
  run_corrected.sh           Driver: validation gate -> main run
  run_lite.R                 Fast preview (fewer chains/iters)
  run_predict_only.R         Prediction-only re-run from a saved fit
  run_prevalence_cv_metrics.R  Spatial block-CV prevalence metrics (AUC, Moran's I)
```

## Data availability

**No data are included in this repository.** The analysis uses individual-level
*Loa loa* survey records (which contain potentially identifying location information)
and gridded environmental covariates. The survey data are available from the study
custodians under the consent and ethics terms of the original survey; the
environmental covariates can be regenerated with `code/NGA_gee_extract.py` (requires
a Google Earth Engine account and your own Cloud project ID).

Scripts expect inputs under `data/` and write to `results/` relative to the
repository root (see `CFG$paths` in `NGA_common.R`):

```
data/NGA_cleaned_corrected.csv                  individual-level survey data (not provided)
data/prediction_NGA_states.csv                  ~2 km prediction grid
data/covariates/covariates_villages_corrected.csv
data/covariates/covariates_grid.csv
```

## Requirements

- **R** (≥ 4.1) with `rstan` (and a working C++ toolchain for Stan), `loo`,
  `tidyverse`/`dplyr`/`readr`/`purrr`, `fields`, `splines`, `pROC`.
- **Python** (≥ 3.9) with `earthengine-api` for covariate extraction (optional).

## Running

Run from the repository root (scripts resolve paths relative to it):

```bash
# 1. validation gate, then the main high-risk run
bash code/run_corrected.sh

# 2. supplementary analyses
Rscript code/NGA_multithreshold.R && Rscript code/NGA_multithreshold_report.R
Rscript code/NGA_binomial_benchmark.R
Rscript code/NGA_report_figures.R

# 3. tail-model adequacy and out-of-sample checks (supplementary material)
Rscript code/NGA_shape_stability.R          # threshold stability, u = 500..1500
Rscript code/NGA_prior_sensitivity_xi.R     # prior sensitivity of the shape xi
Rscript code/NGA_village_cv.R               # random village-fold CV

# 4. spatially-structured extreme-value tail (Peter Diggle's suggestion)
Rscript code/NGA_spatial_tail_prototype.R   # fast two-stage prototype
Rscript code/NGA_joint_tail.R               # full joint fit (set JOINT_SMOKE=1 to test)
Rscript code/NGA_joint_tail_compare.R       # vs the pooled-tail headline

# 5. revised manuscript figures (fast; no refit needed)
Rscript code/NGA_manuscript_figures_rev.R
```

Scripts that produce manuscript figures write to `results/highrisk_20000/figures/`
(or `.../validation/`). Set `LOA_MS_FIGDIR` to mirror them into a manuscript
directory as well:

```bash
LOA_MS_FIGDIR=/path/to/manuscript/figures Rscript code/NGA_manuscript_figures_rev.R
```

Full MCMC settings (4 chains × 2,000 iterations, adapt_delta 0.95) are in `CFG` at
the top of `NGA_common.R`; reduce them (or use `run_lite.R`) for a quick check.

## Input and output overrides

`CFG$paths$data` and `CFG$paths$out_dir` default to
`data/NGA_cleaned_corrected.csv` and `results/highrisk_20000`, and can be pointed
elsewhere with the `LOA_DATA` and `LOA_OUT` environment variables. This allows a
parallel run on a corrected input without disturbing an existing results directory:

```bash
LOA_DATA=data/NGA_lgafix.csv LOA_OUT=results/highrisk_20000_fixed \
  Rscript code/NGA_run_highrisk.R
```

## Village identity and the LGA correction

Villages were originally keyed on the administrative polygon each GPS point falls
in (`NAME_1`, `NAME_2`, `Vill`). Several villages carry points in two different
LGAs — four of them tens of km apart, which are data-entry errors — so those
villages were split in two and the counts (116 villages, 25 LGAs) departed from the
survey's own labelling (110 communities, 19 LGAs). `NGA_lga_fix_build.R` keys on
`VillID` and keeps, for each village, the recorded points whose LGA matches the
survey's own label. `NGA_lga_fix_run.R` fits the headline model on both inputs and
compares the per-LGA results; the safety conclusion is unchanged.

## Requirements note

`build_covariates()` places natural-spline interior knots at quantiles of the
distinct covariate values, dropping any that coincide with a boundary knot;
`water_occurrence` is zero at most villages and a naive quantile rule puts several
interior knots on the left boundary, which `splines::ns()` rejects. Note also that
patchwork must be ≥ 1.3.0 to work with ggplot2 4.x — earlier versions fail when
collecting guides across panels.

## License

MIT — see [LICENSE](LICENSE).
