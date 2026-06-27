#!/usr/bin/env bash
# =============================================================================
# run_corrected.sh — full rerun on the corrected data (Edo->Delta, Benue->Ebonyi).
# Validation gate, then the main high-risk run. projpred/subset are intentionally
# dropped (covariate story is the in-model horseshoe). Run from the Nigeria/ folder.
# =============================================================================
set -u
cd "$(dirname "$0")/.." || exit 1
MASTER=results/corrected_pipeline.log; : > "$MASTER"
stamp() { date +%H:%M; }
run_step() {
  echo "[$(stamp)] START $1" | tee -a "$MASTER"
  if Rscript "$3" > "results/$2" 2>&1; then echo "[$(stamp)] OK    $1" | tee -a "$MASTER"
  else echo "[$(stamp)] FAIL  $1 (see results/$2)" | tee -a "$MASTER"; exit 1; fi
}
run_step validation corrected_1_validation.log code/NGA_validation.R
run_step run        corrected_2_run.log        code/NGA_run_highrisk.R
echo "[$(stamp)] DONE validation + run" | tee -a "$MASTER"
