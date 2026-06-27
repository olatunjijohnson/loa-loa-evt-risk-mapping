# =============================================================================
# run_lite.R — fast APPROXIMATE pass of the full high-risk pipeline.
# Real model + all 25 LGAs + all 4 tails, but reduced MCMC (2 chains, 600 iter,
# adapt_delta 0.85, treedepth 9) and thinned prediction. ~1.5h. Results land in
# results/highrisk_20000/lite/. NOT for final inference — a preview before the
# full rigorous run. Run from the Nigeria/ folder:  Rscript code/run_lite.R
# =============================================================================
RUN_OVERRIDE <- list(
  chains = 2, iter = 500, warmup = 250,
  adapt_delta = 0.85, max_treedepth = 9, thin_pred = 4
)
OUT_DIR_OVERRIDE <- "results/highrisk_20000/lite"
source("code/NGA_run_highrisk.R")
