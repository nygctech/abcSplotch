library(INLA)
library(rstan)
library(Matrix)

source("LCAR.r")

script_start_time <- Sys.time()

########################### READ INPUT DATA ###########################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript model_script.R <input_rdat_path> <output_directory> [TRUE/FALSE (sample posterior)]")
}

input_path  <- args[1]
output_root <- args[2]
draw_samples <- if (length(args) > 2) { tolower(args[4]) %in% c("true", "t", "1", "yes") } else FALSE

# Assumes pattern like: inputs/0/data_1.R
base_name <- basename(input_path)

# Extract number after last underscore
index <- sub(".*_([0-9]+)\\.[Rr].*$", "\\1", base_name)

if (is.na(index) || index == base_name) {
  stop("Could not extract gene index from input filename.")
}

rdat = read_rdump(input_path)

N = sum(rdat$N_spots)
D = rdat$D
y = rdat$counts
size_factors = rdat$size_factors
log_size_factors = log(size_factors)

# Check whether compositional and/or spatial information is present; adjust model accordingly
compositional = ("N_celltypes" %in% names(rdat))
spatial = (rdat$car > 0)

if (compositional) {
    K = rdat$N_celltypes
    E = rdat$E
}

# Create sparse adjacency matrix
if (spatial) {
    i <- rdat$W_sparse[, 1]
    j <- rdat$W_sparse[, 2]
    
    # Build symmetric adjacency directly (still sparse CSC / dgCMatrix)
    W_sparse <- sparseMatrix(
      i = c(i, j),
      j = c(j, i),
      x = 1,
      dims = c(N, N),
      giveCsparse = TRUE
    )
    
    # If duplicates exist, sparseMatrix will sum them (2,3,...) -> clamp back to 1
    W_sparse@x[] <- 1
    
    # Ensure no self edges (optional safety)
    Matrix::diag(W_sparse) <- 0
    W_sparse <- Matrix::drop0(W_sparse)
}

# Derive spot-level labeling for each covariate group
if (rdat$N_level_3 > 0) {
    level_3 = rep(rdat$tissue_mapping, times = rdat$N_spots)
    level_2 = rdat$level_3_mapping[level_3]
    level_1 = rdat$level_2_mapping[level_2]
} else if (rdat$N_level_2 > 0) {
    level_2 = rep(rdat$tissue_mapping, times = rdat$N_spots)
    level_1 = rdat$level_2_mapping[level_2]
} else {
    level_1 = rep(rdat$tissue_mapping, times = rdat$N_spots)
}

########################### CONSTRUCT DESIGN MATRIX ###########################

# Make factors for condition group and MROI
cond_fac = factor(level_1)
mroi_fac = factor(D)

# Combine condition + MROI into a single group index
g_fac = interaction(cond_fac, mroi_fac, drop = TRUE)  # levels correspond to (cond, MROI) combos -- e.g., "1.2"
g_id = as.integer(g_fac)
G = nlevels(g_fac)

if (compositional) {
  stopifnot(nrow(E) == N, ncol(E) == K)

  # For each row i, nonzeros are in columns (g_id[i]-1)*K + 1..K with values E[i,]
  i_idx = rep(seq_len(N), each = K)                        # length N*K
  j_idx = (rep(g_id, each = K) - 1L) * K + rep(seq_len(K), times = N)
  x_val = as.numeric(t(E))                                 # row-wise flatten: E[1,], E[2,], ...

  X = sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(N, G * K),
    giveCsparse = TRUE
  )

} else {
  # One-hot: row i has a single 1 in column g_id[i]
  X = sparseMatrix(
    i = seq_len(N),
    j = g_id,
    x = 1,
    dims = c(N, G),
    giveCsparse = TRUE
  )
}

X = drop0(X)  # optional cleanup

# Make readable coefficient names: beta[g,k] where g encodes (condition, MROI)
g_levels = levels(g_fac)  # by default, strings like "1.2"
parts = strsplit(g_levels, ".", fixed = TRUE)
c_idx = match(vapply(parts, `[`, "", 1), levels(cond_fac))
m_idx = match(vapply(parts, `[`, "", 2), levels(mroi_fac))

stopifnot(length(c_idx) == G, length(m_idx) == G)  # sanity check; c_idx and m_idx should both be equal to number of conditions

if (compositional) {
    colnames(X) <- as.vector(sapply(seq_len(G), function(g) {
        paste0("beta_c", c_idx[g], "_m", m_idx[g], "_k", seq_len(K))
    }))
} else {
    colnames(X) <- as.vector(sapply(seq_len(G), function(g) {
        paste0("beta_c", c_idx[g], "_m", m_idx[g])
    }))
}

# Tell INLA how to map data to parameters
# - inla.stack is a fancy DataFrame that can handle different-dimension input
# - Our "betas" are length p (conds * mrois * celltypes) while spot_ids, scale_factors are length N
p = ncol(X)

effects = list(
    beta = 1:p,
    sample_id = 1:N
)

A = list(
    X,                 # maps beta (length p) -> N observations
    1                  # maps sample_id (length N) -> N observations (identity)
)

if (spatial) {
    effects$region_id = 1:N
    A = c(A, list(1))       # maps region_id (length N) -> N observations (identity)
}

stk = inla.stack(
    data = list(y = y),
    A = A,
    effects = effects,
    tag = "est"
)

########################### DEFINE PRIORS ###########################

# BETA:
beta_prec = 0.25   # SD = 2, as in Stan formulation
# TODO: allow modulation of beta_mean and beta_prec?

# EPSILON:
# In Stan implementation, hyperprior on SD of epsilon is HalfNormal(0,0.3)
# These PC prior settings make sure that SD of epsilon rarely exceeds what we'd expect from this
hyper_iid <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.773, 0.01)  # P(sigma > 0.773) = 0.01
  )
)

# PSI: (within LCAR.r)

########################### BUILD & RUN MODEL ###########################

# Formula: no intercept; offset passed later to fit; beta is a latent coefficient vector
fml = y ~ 0 +
  f(beta, model = "iid",
    hyper = list(
      prec = list(initial = log(beta_prec), fixed = TRUE)
    )
  ) +
  f(sample_id, model = "iid", hyper = hyper_iid)

if (spatial) {
    # Define Leroux CAR component
    LCAR.model = inla.LCAR.model(W=W_sparse)
    fml = update(fml, . ~ . + f(region_id, model = LCAR.model))
}

fit_args = list(
    formula = fml,
    family = "nbinomial",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    offset = log_size_factors
)
fit = do.call(INLA::inla, fit_args)

# Quick check
beta_summary = fit$summary.random$beta
beta_summary$coef = colnames(X)
print(beta_summary[1:min(10, nrow(beta_summary)), ])
print(fit$summary.hyperpar)

########################### SAMPLE JOINT POSTERIOR ###########################

# Expensive; and not yielding comparable results to BF testing on Bayesian posteriors.
if (draw_samples) {
    # Sample joint posterior, only saving results from betas:
    cat("Sampling joint posterior from computed marginals...\n")
    
    # Chunking reduces peak memory usage, as INLA generates full posterior then returns subset
    n_total <- 1000
    chunk_size <- 50
    
    beta_samples <- matrix(NA_real_, n_total, p)
    
    starts <- seq(1, n_total, by = chunk_size)
    
    for (s in starts) {
      m <- min(chunk_size, n_total - s + 1)
    
      samp_chunk <- inla.posterior.sample(
        n = m,
        result = fit,
        selection = list(beta = 1:p),
        add.names = FALSE
      )
    
      beta_samples[s:(s + m - 1), ] <- t(vapply(
        samp_chunk,
        function(z) as.numeric(z$latent[, 1]),
        numeric(p)
      ))
    
      rm(samp_chunk)
      gc()
    }
    
    colnames(beta_samples) <- colnames(X)
}

########################### SAVE RESULTS ###########################

# Full plot-ready marginals for beta (now in marginals.random$beta)
beta_marginals <- fit$marginals.random$beta
names(beta_marginals) <- colnames(X)

# Separate per-spot/cell estimates (defined as "est" in inla.stack) from fixed effects (betas)
idx_est = inla.stack.index(stk, tag = "est")$data
stopifnot(length(idx_est) == length(size_factors))

# Full plot-ready marginals & joint posteriors for Beta; summary statistics for other marginals & hyperparameters
res = list(
  # store the controls used (so later scripts know what was requested)
  control = list(
    predictor = list(compute = TRUE),
    compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  ),

  # summaries
  beta_summary  = beta_summary,
  hyper_summary  = fit$summary.hyperpar,

  log_lambda = list(
      mean = fit$summary.linear.predictor$mean[idx_est],
      sd = fit$summary.linear.predictor$sd[idx_est]
  ),
  lambda = list(
      mean = fit$summary.fitted.values$mean[idx_est] / size_factors,
      sd = fit$summary.fitted.values$sd[idx_est] / size_factors
  ),
  epsilon = list(
      mean = fit$summary.random$sample_id$mean,
      sd = fit$summary.random$sample_id$sd
  ),

  # keep only beta marginals to keep file size down
  beta_marginals = beta_marginals
)
if (spatial) {
    res$psi = list(
        mean = fit$summary.random$region_id$mean,
        sd = fit$summary.random$region_id$sd
    )
}
if (draw_samples) {
    res$beta_samples = beta_samples
}

# Summaries of model evaluation criteria: DIC, WAIC, CPO
res$criteria = list(
  dic = list(
    value  = fit$dic$dic,        # DIC score (lower is better)
    p_eff  = fit$dic$p.eff       # Effective number of parameters (model complexity)
  ),
  waic = list(
    value  = fit$waic$waic,      # WAIC score (preferred over DIC; lower is better)
    p_eff  = fit$waic$p.eff      # Effective number of parameters under WAIC
  ),
  cpo = list(
    sum_lcpo = sum(fit$cpo$lcpo, na.rm = TRUE),   # Sum log CPO (higher is better predictive fit)
    n_fail   = sum(fit$cpo$failure)               # Number of failed CPO evaluations (should be 0)
  )
)

# ---- Create output directory ----
subdir = as.character(floor(as.integer(index)/100))
output_dir <- file.path(output_root, subdir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Construct output filename ----
output_file <- file.path(output_dir, paste0("results_", index, ".R"))

saveRDS(res, file = output_file, compress = "xz")

########################### RUNTIME LOGGING ###########################

cat("Saved results to:", output_file, "\n")
cat("Wrote file size (bytes):", file.info(output_file)$size, "\n")

script_end_time <- Sys.time()
elapsed_secs <- as.numeric(difftime(script_end_time, script_start_time, units = "secs"))
hours <- floor(elapsed_secs / 3600)
minutes <- floor((elapsed_secs %% 3600) / 60)
seconds <- elapsed_secs %% 60

cat("\n")
cat("Script started:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Script ended:  ", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("Runtime: %02d:%02d:%05.2f (hh:mm:ss)\n", hours, minutes, seconds))