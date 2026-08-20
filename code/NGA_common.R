# =============================================================================
# NGA_common.R
# Shared setup for the Nigeria Loa loa "high-risk (>20,000 mf/mL)" analysis.
#   - Data loading / cleaning (5 target states)
#   - Village-level spatial index + distance matrix
#   - Covariate design (placeholder until GEE covariates land)
#   - Stan model code: JOINT coregionalized GP hurdle (Weibull / lognormal / gamma)
#                      + pooled generalized-Pareto (POT/EVT) tail
#   - Helper functions: fit, krige to 2 km grid, exceedance, LGA aggregation
#
# Sourced by NGA_run_highrisk.R (answer) and NGA_validation.R (gating checks).
# Coordinates are UTM in KILOMETRES (utmx, utmy) -> GP range rho is in km.
# NOTE: nothing here is run on source except function/string definitions.
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(rstan)
})
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
CFG <- list(
  target_states = c("Abia", "Anambra", "Delta", "Ebonyi", "Imo"),
  hr_threshold  = 20000,     # mf/mL "high-risk" threshold of interest
  ix_threshold  = 0.001,     # 0.1% — collaborator's "<0.1% of individuals"
  gpd_u         = 1000,      # POT threshold for the EVT tail (revisited in validation)
  # MCMC (full run, feasible-but-rigorous). lite/CV/test override these.
  chains        = 4,
  iter          = 2000,
  warmup        = 1000,
  adapt_delta   = 0.95,
  max_treedepth = 12,
  seed          = 1234,
  thin_pred     = 2,         # thin posterior before exact per-LGA prediction
  village_key   = c("NAME_1", "NAME_2", "Vill"),  # defines a "village" location
  use_covariates = TRUE,     # FALSE -> intercept-only (no GEE covariates)
  cov_subset    = NULL,      # NULL -> all COVARS; else a subset (e.g. projpred pick)
  spline_df     = 4,         # natural-spline df per continuous covariate
  hs_p0         = 6,         # regularized horseshoe: prior guess at # relevant terms
  paths = list(
    # Defaults unchanged; LOA_DATA / LOA_OUT allow a parallel run on a corrected
    # input without disturbing results/highrisk_20000 (see run_lgafix_full.sh).
    data      = Sys.getenv("LOA_DATA", "data/NGA_cleaned_corrected.csv"),  # Edo->Delta, Benue->Ebonyi (see NGA_fix_data.R)
    grid      = "data/prediction_NGA_states.csv",
    cov_vill  = "data/covariates/covariates_villages_corrected.csv",  # +2 villages, nearest-cell covariates
    cov_grid  = "data/covariates/covariates_grid.csv",
    out_dir   = Sys.getenv("LOA_OUT", "results/highrisk_20000")
  )
)

# Continuous environmental covariates (from Earth Engine, 2016).
COVARS <- c("elevation", "slope", "lst_day", "ndvi", "rain_total",
            "rain_seasonality", "treecover", "water_occurrence", "population")

# Map family name -> integer flag used inside the Stan hurdle model.
FAMILY_CODE <- c(weibull = 1L, lognormal = 2L, gamma = 3L)

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
# Returns a list with everything the models need:
#   data_model : individual rows (5 states), with y_binary / y_intensity / vid
#   villages   : one row per unique village location (key + mean utmx/utmy)
#   Dmat       : V x V village distance matrix (km)
#   grid       : prediction grid (2 km) restricted to the 5 states
#   Dgv        : N_grid x V cross-distance matrix (km) for kriging
#   X, X_grid  : covariate design matrices (0 columns until GEE covariates added)
load_nga <- function(cfg = CFG) {
  d <- readr::read_csv(cfg$paths$data, show_col_types = FALSE) |>
    filter(NAME_1 %in% cfg$target_states) |>
    filter(!is.na(utmx), !is.na(utmy), !is.na(Loa_count)) |>
    mutate(
      y_binary    = as.integer(Loa_count > 0),
      y_intensity = as.numeric(Loa_count)
    )

  # Village = unique location key; its coordinate is the mean of member rows.
  d <- d |>
    unite("vill_key", all_of(cfg$village_key), remove = FALSE, sep = " | ")
  villages <- d |>
    group_by(vill_key) |>
    summarise(utmx = mean(utmx), utmy = mean(utmy),
              NAME_1 = first(NAME_1), NAME_2 = first(NAME_2),
              Vill = first(Vill), .groups = "drop") |>
    arrange(vill_key) |>
    mutate(vid = row_number())
  d <- d |> left_join(villages |> select(vill_key, vid), by = "vill_key")

  V <- nrow(villages)
  vcoords <- as.matrix(villages[, c("utmx", "utmy")])
  Dmat <- as.matrix(dist(vcoords))                       # V x V (km)

  use_cov <- cfg$use_covariates && file.exists(cfg$paths$cov_grid)

  # Prediction grid (2 km, 5 states). With covariates we read the covariate
  # grid file directly so grid rows and covariate rows are inherently aligned.
  if (use_cov) {
    grid <- readr::read_csv(cfg$paths$cov_grid, show_col_types = FALSE) |>
      filter(NAME_1 %in% cfg$target_states)
  } else {
    grid <- readr::read_csv(cfg$paths$grid, show_col_types = FALSE) |>
      filter(NAME_1 %in% cfg$target_states)
  }
  gcoords <- as.matrix(grid[, c("utmx", "utmy")])
  Dgv <- fields_rdist(gcoords, vcoords)                  # N_grid x V (km)

  # Covariate design: natural-spline basis (non-linearity), village-level,
  # mapped to individuals via vid and to the grid cells.
  if (use_cov) {
    cb <- build_covariates(villages, grid, cfg)
    Xv <- cb$Xv; rownames(Xv) <- villages$vill_key   # village-level basis (for CV)
    X      <- Xv[d$vid, , drop = FALSE]              # individual = his village's basis
    X_grid <- cb$Xg
    cov_meta <- cb$meta
  } else {
    Xv <- matrix(0, V, 0); X <- matrix(0, nrow(d), 0)
    X_grid <- matrix(0, nrow(grid), 0); cov_meta <- NULL
  }

  list(d = d, villages = villages, V = V, Dmat = Dmat,
       grid = grid, Dgv = Dgv, X = X, X_grid = X_grid, Xv = Xv, cov_meta = cov_meta,
       pos_idx = which(d$y_binary == 1))
}

# -----------------------------------------------------------------------------
# Build covariate spline-basis design matrices (villages + grid).
#   - population log1p-transformed; each covariate standardized on village stats
#   - natural cubic spline (ns, df) per covariate; same knots applied to grid
#   - returns Xv (V x P), Xg (N_grid x P), and column->covariate metadata
# -----------------------------------------------------------------------------
build_covariates <- function(villages, grid, cfg = CFG) {
  covs <- if (!is.null(cfg$cov_subset)) cfg$cov_subset else COVARS  # full set or subset
  cvill <- readr::read_csv(cfg$paths$cov_vill, show_col_types = FALSE)
  prep  <- function(D) dplyr::mutate(D,
            population = log1p(ifelse(is.na(population), 0, population)))
  cvill <- prep(cvill)
  vjoin <- villages |>
    left_join(cvill |> select(NAME_1, NAME_2, Vill, all_of(covs)),
              by = c("NAME_1", "NAME_2", "Vill"))
  if (any(is.na(vjoin[covs])))
    stop("Missing covariate values for some villages after join.")
  gridc <- prep(grid)

  mu <- vapply(covs, function(c) mean(vjoin[[c]]), numeric(1))
  sg <- vapply(covs, function(c) sd(vjoin[[c]]),   numeric(1))

  Bv <- list(); Bg <- list(); meta <- character()
  for (c in covs) {
    xv <- (vjoin[[c]] - mu[c]) / sg[c]
    xg <- (gridc[[c]] - mu[c]) / sg[c]
    # Interior knots at quantiles of the DISTINCT values, dropping any that fall on
    # a boundary knot. water_occurrence is zero at most villages, so quantiles of
    # all values put several interior knots on the left boundary, which ns() rejects.
    bk <- range(xv)
    k  <- unique(as.numeric(quantile(unique(xv),
            seq(0, 1, length.out = cfg$spline_df + 1)[-c(1, cfg$spline_df + 1)])))
    k  <- k[k > bk[1] & k < bk[2]]
    if (length(k)) {
      b <- splines::ns(xv, knots = k, Boundary.knots = bk)
      Bv[[c]] <- b
      Bg[[c]] <- predict(b, xg)
    } else {                                   # degenerate: fall back to linear
      Bv[[c]] <- matrix(xv, ncol = 1)
      Bg[[c]] <- matrix(xg, ncol = 1)
      message("  build_covariates: linear fallback for ", c)
    }
    meta <- c(meta, rep(c, ncol(Bv[[c]])))
  }
  Xv <- do.call(cbind, Bv); Xg <- do.call(cbind, Bg)
  colnames(Xv) <- colnames(Xg) <- paste0(meta, ".", ave(meta, meta, FUN = seq_along))
  list(Xv = Xv, Xg = Xg, meta = meta, mu = mu, sd = sg)
}

# Lightweight Euclidean cross-distance (avoids a hard dependency on fields::rdist)
fields_rdist <- function(A, B) {
  # ||a-b||^2 = |a|^2 + |b|^2 - 2 a.b
  aa <- rowSums(A^2); bb <- rowSums(B^2)
  d2 <- outer(aa, bb, "+") - 2 * (A %*% t(B))
  sqrt(pmax(d2, 0))
}

# -----------------------------------------------------------------------------
# Stan model: JOINT GP hurdle (Weibull / lognormal / gamma) with a BIVARIATE
# intrinsic coregionalization (ICM) on the two spatial fields.
#   prevalence:  logit pi  = b0_p + X beta_p + f_p(s)
#   intensity:   link(.)   = b0_c + X beta_c + f_c(s)   (positives only)
#   spatial:     (f_p, f_c)(s) ~ GP with separable covariance Sigma (x) K(rho),
#                K(rho) = exp(-d/rho) a single shared spatial correlation, and
#                Sigma = diag(sigma_f) Omega diag(sigma_f) the 2x2 cross-covariance.
#   correlation(f_p, f_c) = Omega[1,2] (reported as rho_pc).
#   NB: this REPLACES the earlier shared-factor + loadings form (w, lambda_p,
#   lambda_c). There both loadings were ~0 and unidentified -> a funnel that gave
#   R-hat ~1.14 / ESS ~30 / divergences. The ICM parameterizes the cross-corr
#   directly (LKJ prior) and is well-identified.
# -----------------------------------------------------------------------------
STAN_HURDLE <- "
data {
  int<lower=1> N;
  int<lower=1> V;
  array[N] int<lower=1,upper=V> vid;
  array[N] int<lower=0,upper=1> y;
  int<lower=0> Npos;
  array[Npos] int<lower=1,upper=N> pos_idx;
  vector<lower=0>[Npos] yc;
  int<lower=0> P;
  matrix[N, P] X;
  matrix[V, V] Dmat;
  int<lower=1,upper=3> family;   // 1 weibull, 2 lognormal, 3 gamma
  real<lower=0> hs_scale_global; // regularized-horseshoe global scale (tau0)
}
transformed data {
  real jitter = 1e-8;
  real hs_slab_scale = 2;        // slab scale c
  real hs_slab_df = 4;           // slab degrees of freedom
}
parameters {
  real b0_p;
  real b0_c;

  // regularized horseshoe — prevalence (covariate) coefficients
  vector[P] zb_p;
  vector<lower=0>[P] hs_lam_p;
  real<lower=0> hs_tau_p;
  real<lower=0> hs_caux_p;
  // regularized horseshoe — intensity (covariate) coefficients
  vector[P] zb_c;
  vector<lower=0>[P] hs_lam_c;
  real<lower=0> hs_tau_c;
  real<lower=0> hs_caux_c;

  // bivariate intrinsic coregionalization (ICM): one shared spatial range and a
  // 2x2 cross-covariance Sigma = diag(sigma_f) Omega diag(sigma_f).
  real<lower=0> rho;                 // shared spatial range (km)
  vector<lower=0>[2] sigma_f;        // marginal field SDs: [1]=prevalence, [2]=intensity
  cholesky_factor_corr[2] L_Omega;   // correlation Cholesky (LKJ prior)
  matrix[V, 2] Z;                    // std-normal spatial innovations

  real<lower=0> shape;     // weibull alpha OR gamma shape
  real<lower=0> sigma_ln;  // lognormal sd
}
transformed parameters {
  vector[V] f_p;
  vector[V] f_c;
  vector[P] beta_p;
  vector[P] beta_c;
  {
    real c2p = square(hs_slab_scale) * hs_caux_p;
    real c2c = square(hs_slab_scale) * hs_caux_c;
    vector[P] lt_p = sqrt(c2p * square(hs_lam_p) ./ (c2p + square(hs_tau_p) * square(hs_lam_p)));
    vector[P] lt_c = sqrt(c2c * square(hs_lam_c) ./ (c2c + square(hs_tau_c) * square(hs_lam_c)));
    beta_p = zb_p .* lt_p * hs_tau_p;
    beta_c = zb_c .* lt_c * hs_tau_c;
  }
  {
    // ICM: vec((f_p, f_c)) has covariance Sigma (x) K(rho).  F = Lk Z L_Sigma'
    // gives Cov = (L_Sigma L_Sigma') (x) (Lk Lk') = Sigma (x) K.
    matrix[V, V] Lk = cholesky_decompose(add_diag(exp(-Dmat / rho), jitter));
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_f, L_Omega);  // Cholesky of Sigma
    matrix[V, 2] F = Lk * Z * L_Sigma';
    f_p = F[, 1];
    f_c = F[, 2];
  }
}
model {
  // priors
  b0_p ~ normal(0, 2);
  b0_c ~ normal(0, 2);
  // regularized horseshoe priors (Piironen & Vehtari)
  zb_p ~ std_normal();
  zb_c ~ std_normal();
  hs_lam_p ~ cauchy(0, 1);
  hs_lam_c ~ cauchy(0, 1);
  hs_tau_p ~ cauchy(0, hs_scale_global);
  hs_tau_c ~ cauchy(0, hs_scale_global);
  hs_caux_p ~ inv_gamma(0.5 * hs_slab_df, 0.5 * hs_slab_df);
  hs_caux_c ~ inv_gamma(0.5 * hs_slab_df, 0.5 * hs_slab_df);
  rho ~ inv_gamma(5, 150);          // range prior centred ~ 30-40 km (km units)
  sigma_f ~ normal(0, 1);           // half-normal marginal field SDs
  L_Omega ~ lkj_corr_cholesky(2);   // mild shrinkage of the prev-intensity corr toward 0
  to_vector(Z) ~ std_normal();
  shape ~ normal(1, 0.5);
  sigma_ln ~ normal(0, 1);

  // prevalence
  {
    vector[N] eta_p = b0_p + f_p[vid];
    if (P > 0) eta_p += X * beta_p;
    y ~ bernoulli_logit(eta_p);
  }

  // intensity (positives only)
  for (k in 1:Npos) {
    int i = pos_idx[k];
    real eta = b0_c + f_c[vid[i]];
    if (P > 0) eta += X[i] * beta_c;
    if (family == 1)      yc[k] ~ weibull(shape, exp(eta));
    else if (family == 2) yc[k] ~ lognormal(eta, sigma_ln);
    else                  yc[k] ~ gamma(shape, shape / exp(eta));
  }
}
generated quantities {
  real rho_pc = multiply_lower_tri_self_transpose(L_Omega)[2, 1];  // prev-intensity spatial corr
  vector[N] log_lik_bin;
  vector[Npos] log_lik_int;
  {
    vector[N] eta_p = b0_p + f_p[vid];
    if (P > 0) eta_p += X * beta_p;
    for (i in 1:N) log_lik_bin[i] = bernoulli_logit_lpmf(y[i] | eta_p[i]);
  }
  for (k in 1:Npos) {
    int i = pos_idx[k];
    real eta = b0_c + f_c[vid[i]];
    if (P > 0) eta += X[i] * beta_c;
    if (family == 1)      log_lik_int[k] = weibull_lpdf(yc[k] | shape, exp(eta));
    else if (family == 2) log_lik_int[k] = lognormal_lpdf(yc[k] | eta, sigma_ln);
    else                  log_lik_int[k] = gamma_lpdf(yc[k] | shape, shape / exp(eta));
  }
}
"

# -----------------------------------------------------------------------------
# Stan model: pooled generalized-Pareto tail (peaks-over-threshold)
#   For exceedances e = (count - u) | count > u :  e ~ GPD(sigma, xi)
#   Survival ratio used downstream:
#     R(T) = P(count > T | count > u) = (1 + xi (T-u)/sigma)^(-1/xi)   (xi != 0)
# This estimates the TAIL SHAPE only; the spatially-varying P(count>u | x)
# is taken from the spatial body model (see exceedance_gpd()).
# -----------------------------------------------------------------------------
STAN_GPD <- "
functions {
  real gpd_lpdf(vector e, real sigma, real xi) {
    int n = num_elements(e);
    if (fabs(xi) < 1e-8)
      return -n * log(sigma) - sum(e) / sigma;        // xi -> 0 : exponential
    else {
      vector[n] z = 1 + xi * e / sigma;
      for (i in 1:n) if (z[i] <= 0) return negative_infinity();
      return -n * log(sigma) - (1 + 1 / xi) * sum(log(z));
    }
  }
}
data {
  int<lower=1> M;            // number of exceedances over u
  vector<lower=0>[M] e;      // count - u
}
parameters {
  real<lower=0> sigma;
  real<lower=-0.5> xi;       // restrict to finite-variance-ish tails
}
model {
  sigma ~ normal(0, 2000);
  xi ~ normal(0, 0.5);
  e ~ gpd(sigma, xi);
}
"

# -----------------------------------------------------------------------------
# Fitting helpers
# -----------------------------------------------------------------------------
compile_models <- function() {
  list(
    hurdle = stan_model(model_code = STAN_HURDLE),
    gpd    = stan_model(model_code = STAN_GPD)
  )
}

fit_hurdle <- function(dat, family = c("weibull", "lognormal", "gamma"),
                       model, cfg = CFG, iter = cfg$iter, warmup = cfg$warmup,
                       chains = cfg$chains) {
  family <- match.arg(family)
  P  <- ncol(dat$X)
  N  <- nrow(dat$d)
  p0 <- cfg$hs_p0
  tau0 <- if (P > p0) p0 / (P - p0) / sqrt(N) else 1  # horseshoe global scale
  stan_data <- list(
    N = N, V = dat$V, vid = dat$d$vid, y = dat$d$y_binary,
    Npos = length(dat$pos_idx), pos_idx = dat$pos_idx,
    yc = dat$d$y_intensity[dat$pos_idx],
    P = P, X = dat$X, Dmat = dat$Dmat,
    family = FAMILY_CODE[[family]], hs_scale_global = tau0
  )
  sampling(model, data = stan_data, chains = chains, iter = iter, warmup = warmup,
           seed = cfg$seed, control = list(adapt_delta = cfg$adapt_delta,
                                           max_treedepth = cfg$max_treedepth))
}

fit_gpd <- function(dat, model, cfg = CFG) {
  yc <- dat$d$y_intensity[dat$pos_idx]
  e  <- yc[yc > cfg$gpd_u] - cfg$gpd_u
  sampling(model, data = list(M = length(e), e = e),
           chains = cfg$chains, iter = cfg$iter, warmup = cfg$warmup, seed = cfg$seed)
}

# -----------------------------------------------------------------------------
# Per-cell tail survival P(count > T | infected): scalar shape/sd, vector eta.
# -----------------------------------------------------------------------------
surv_cells <- function(family, eta, shape = NULL, sigma_ln = NULL, Tval) {
  switch(family,
    weibull   = exp(-(Tval / exp(eta))^shape),
    lognormal = plnorm(Tval, meanlog = eta, sdlog = sigma_ln, lower.tail = FALSE),
    gamma     = pgamma(Tval, shape = shape, rate = shape / exp(eta), lower.tail = FALSE)
  )
}

# GPD conditional survival ratio R(T) = P(count > T | count > u) for one draw.
gpd_ratio <- function(sigma, xi, u, Tval) {
  if (abs(xi) < 1e-8) exp(-(Tval - u) / sigma)
  else max(1 + xi * (Tval - u) / sigma, 0)^(-1 / xi)
}

# -----------------------------------------------------------------------------
# Mean kriging (conditional MEAN only) — used for CV holdout prediction & maps.
#   z_draws: S x V latent village field; returns S x Nx conditional means.
# -----------------------------------------------------------------------------
krige_mean <- function(z_draws, Dmat, Dxv, rho_draws, sig_draws = NULL) {
  S <- nrow(z_draws); Nx <- nrow(Dxv); out <- matrix(0, S, Nx)
  for (s in seq_len(S)) {
    rho <- rho_draws[s]; sg2 <- if (is.null(sig_draws)) 1 else sig_draws[s]^2
    Kvv <- sg2 * exp(-Dmat / rho); diag(Kvv) <- diag(Kvv) + 1e-8
    Kxv <- sg2 * exp(-Dxv / rho)
    out[s, ] <- as.numeric(Kxv %*% solve(Kvv, z_draws[s, ]))
  }
  out
}
krige_field <- krige_mean   # backward-compatible alias

# -----------------------------------------------------------------------------
# EXACT coregionalized per-LGA prediction (user-requested "exact draw").
#   For each posterior draw and each LGA we draw the bivariate field (f_p, f_c)
#   JOINTLY over that LGA's grid cells from its ICM conditional distribution: the
#   shared spatial correlation K(rho) gives common kriging weights for both
#   components, and the 2x2 Sigma (via L_Sigma) carries the cross-correlation.
#   q_g = mean_cells [ pi(x) * P(count > T | infected, x) ].
#   Per-LGA blocking is EXACT for each LGA's marginal q_g distribution (the
#   cross-LGA field covariance does not affect a within-LGA average).
# -----------------------------------------------------------------------------
precompute_lga_geometry <- function(dat) {
  lga    <- interaction(dat$grid$NAME_1, dat$grid$NAME_2, drop = TRUE)
  vc     <- as.matrix(dat$villages[, c("utmx", "utmy")])
  groups <- split(seq_len(nrow(dat$grid)), lga)
  map(groups, function(cells) {
    gc <- as.matrix(dat$grid[cells, c("utmx", "utmy")])
    # population weights for the per-LGA mean (the answer is about people, not area).
    w  <- dat$grid$population[cells]; w[is.na(w) | w < 0] <- 0
    if (sum(w) <= 0) w <- rep(1, length(cells))   # fall back to area mean if no pop
    list(cells = cells, pop = w,
         NAME_1 = dat$grid$NAME_1[cells[1]], NAME_2 = dat$grid$NAME_2[cells[1]],
         Dgv = fields_rdist(gc, vc), Dgg = as.matrix(dist(gc)))
  })
}

# One JOINT conditional draw of the bivariate field (f_p, f_c) over an LGA's
# cells under the ICM. Spatial correlation K(rho) is shared by both fields, so
# the kriging weights A are common; the 2x2 cross-structure enters via L_Sigma
# (lower Cholesky of Sigma). Returns list(fp, fc), each of length n_g.
.cond_draw_icm <- function(fp_v, fc_v, Kvv_inv, Dgv, Dgg, rho, L_Sigma) {
  Kgv <- exp(-Dgv / rho)                  # unit-variance correlation (Sigma carries variance)
  A   <- Kgv %*% Kvv_inv                  # n_g x V kriging weights (same for both fields)
  mp  <- as.numeric(A %*% fp_v)           # conditional means
  mc  <- as.numeric(A %*% fc_v)
  C   <- exp(-Dgg / rho) - A %*% t(Kgv)   # residual spatial correlation (n_g x n_g)
  diag(C) <- diag(C) + 1e-7
  Lc  <- tryCatch(t(chol(C)), error = function(e) t(chol(C + diag(1e-4, nrow(C)))))
  E   <- matrix(rnorm(nrow(C) * 2), nrow(C), 2)
  resid <- Lc %*% E %*% t(L_Sigma)        # Cov(vec resid) = Sigma (x) C
  list(fp = mp + resid[, 1], fc = mc + resid[, 2])
}

# Returns an S x n_LGA matrix of q_g draws (LGA exceedance proportions).
# tail in {weibull, lognormal, gamma, gpd}; gpd uses the lognormal body for
# P(count>u | x) then applies the pooled EVT ratio R(T).
predict_lga_exact <- function(post, dat, geom, tail, Tval, cfg = CFG,
                              gpd_post = NULL, thin = cfg$thin_pred, seed = cfg$seed) {
  set.seed(seed)
  S0  <- length(post$b0_p); idx <- seq(1, S0, by = thin); S <- length(idx)
  Dmat <- dat$Dmat
  Pc   <- ncol(dat$X_grid)
  body <- if (tail == "gpd") "lognormal" else tail
  qg   <- matrix(0, S, length(geom))
  for (ii in seq_along(idx)) {
    s <- idx[ii]
    rho <- post$rho[s]
    Kvv <- exp(-Dmat / rho); diag(Kvv) <- diag(Kvv) + 1e-8
    Kvv_i <- solve(Kvv)
    # L_Sigma = diag(sigma_f) %*% L_Omega = lower Cholesky of Sigma (matches Stan).
    L_Sigma <- diag(c(post$sigma_f[s, 1], post$sigma_f[s, 2])) %*% post$L_Omega[s, , ]
    R <- if (tail == "gpd") gpd_ratio(gpd_post$sigma[s], gpd_post$xi[s], cfg$gpd_u, Tval) else NA
    for (g in seq_along(geom)) {
      G    <- geom[[g]]
      cd    <- .cond_draw_icm(post$f_p[s, ], post$f_c[s, ], Kvv_i, G$Dgv, G$Dgg, rho, L_Sigma)
      eta_p <- post$b0_p[s] + cd$fp
      eta_c <- post$b0_c[s] + cd$fc
      if (Pc > 0) {
        Xg    <- dat$X_grid[G$cells, , drop = FALSE]
        eta_p <- eta_p + as.numeric(Xg %*% post$beta_p[s, ])
        eta_c <- eta_c + as.numeric(Xg %*% post$beta_c[s, ])
      }
      pi_g  <- plogis(eta_p)
      if (tail == "gpd") {
        su <- surv_cells("lognormal", eta_c, sigma_ln = post$sigma_ln[s], Tval = cfg$gpd_u)
        ix <- pi_g * su * R
      } else {
        sv <- surv_cells(body, eta_c, shape = post$shape[s],
                         sigma_ln = post$sigma_ln[s], Tval = Tval)
        ix <- pi_g * sv
      }
      qg[ii, g] <- weighted.mean(ix, G$pop)   # population-weighted LGA mean
    }
  }
  qg
}

# Summarise q_g draws -> the per-LGA answer: P(q_g < 0.1%) + posterior summary.
summarise_qg <- function(qg, geom, cfg = CFG) {
  map_dfr(seq_along(geom), function(g) {
    q <- qg[, g]
    tibble(
      NAME_1 = geom[[g]]$NAME_1, NAME_2 = geom[[g]]$NAME_2,
      n_cells = length(geom[[g]]$cells),
      q_mean = mean(q), q_lo = quantile(q, .025), q_hi = quantile(q, .975),
      p_below_0.1pct = mean(q < cfg$ix_threshold)   # <-- the collaborator's answer
    )
  })
}
