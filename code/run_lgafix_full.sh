#!/usr/bin/env bash
# Full headline pipeline on the LGA-corrected survey file (data/NGA_lgafix.csv).
# Writes results/highrisk_20000_fixed/; results/highrisk_20000/ is left untouched.
set -euo pipefail
cd "$(dirname "$0")/.."
export LOA_DATA=data/NGA_lgafix.csv
export LOA_OUT=results/highrisk_20000_fixed
mkdir -p "$LOA_OUT/figures"   # multithreshold_report writes here before report_figures creates it
echo "### 1/4 main run  $(date)"; Rscript code/NGA_run_highrisk.R
echo "### 2/4 multithreshold  $(date)"; Rscript code/NGA_multithreshold.R
echo "### 3/4 multithreshold report  $(date)"; Rscript code/NGA_multithreshold_report.R
echo "### 4/4 report figures  $(date)"; Rscript code/NGA_report_figures.R
echo "### complete  $(date)"
