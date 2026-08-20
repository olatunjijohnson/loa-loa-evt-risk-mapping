#!/usr/bin/env bash
# Validation + village-level CV on the LGA-corrected data.
set -uo pipefail
cd "$(dirname "$0")/.."   # repository root
export LOA_DATA=data/NGA_lgafix.csv
export LOA_OUT=results/highrisk_20000_fixed
mkdir -p "$LOA_OUT/validation" "$LOA_OUT/figures"
echo "### validation start $(date)"
Rscript code/NGA_validation.R; echo "### validation exit=$? $(date)"
echo "### village CV start $(date)"
Rscript code/NGA_village_cv.R;  echo "### village CV exit=$? $(date)"
echo "### all validation complete $(date)"
