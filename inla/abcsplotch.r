library(INLA)
library(rstan)
library(Matrix)

source("LCAR.r")

########################### READ INPUT DATA ###########################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript model_script.R <input_rdat_path> <output_directory>")
}

input_path  <- args[1]
output_root <- args[2]

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
    i = rdat$W_sparse[, 1]
    j = rdat$W_sparse[, 2]
    W_sparse = sparseMatrix(
        i = i,
        j = j,
        x = 1,
        dims = c(N, N)
    )
    # sanity enforcement of symmetry
    W_sparse = (W_sparse + t(W_sparse)) > 0
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
    # Build fixed-effect design matrix X of size N x (G*K)
    X = matrix(0, nrow = N, ncol = G * K)
    # ...then populate with E (e.g., for sample i, X[i, (g_i, k)] = E[i,k]
    for (i in seq_len(N)) {
        idx = (g_id[i] - 1) * K + seq_len(K)
        X[i, idx] = E[i,]
    }
} else {
    # Build fixed-effect design matrix X of size N x G
    X = matrix(0, nrow = N, ncol = G)
    for (i in seq_len(N)) {
        X[i, g_id[i]] = 1
    }
}

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

# Dataframe for INLA model
dat = data.frame(
    y = y,
    sample_id = seq_len(N),
    log_size_factor = log_size_factors
)
if (spatial) {
    dat$region_id = seq_len(N)
}

# Bind X into Dataframe (INLA will treat columns as regular covariates)
dat = cbind(dat, as.data.frame(X))

########################### BUILD & RUN MODEL ###########################

if (spatial) {
    # Define Leroux CAR component
    LCAR.model = inla.LCAR.model(W=W_sparse)
    
    # Build GLM
    beta_terms = paste(colnames(X), collapse = " + ")
    glm = as.formula(paste0("y ~ 0 + ", beta_terms, 
                            " + f(sample_id, model='iid', hyper=hyper_iid)",
                            " + f(region_id, model=LCAR.model)",
                            " + offset(log_size_factor)"))
} else {
    # Build GLM
    beta_terms = paste(colnames(X), collapse = " + ")
    glm = as.formula(paste0("y ~ 0 + ", beta_terms, 
                            " + f(sample_id, model='iid', hyper=hyper_iid)",
                            " + offset(log_size_factor)"))
}

# Prior definitions
# BETA:
control.fixed = list(
    mean = list(default = 0), 
    prec = list(default = 0.25)  # SD=2 in Splotch
)
# can be replaced with beta_prior_mean and beta_prior_std vectors of length ncol(X) when snRNA-seq priors available

# EPSILON:
# In Stan implementation, hyperprior on SD of epsilon is HalfNormal(0,0.3)
# These PC prior settings make sure that SD of epsilon rarely exceeds what we'd expect from this
hyper_iid <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.773, 0.01)  # P(sigma > 0.773) = 0.01
  )
)

fit = inla(
    glm,
    family = "nbinomial",
    data = dat,
    control.predictor = list(compute = TRUE),  # save the lambdas (actually nu = log_lambda + log_size_factor)
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.fixed = control.fixed
)

# Quick check
print(fit$summary.fixed[1:min(10, nrow(fit$summary.fixed)), ])
print(fit$summary.hyperpar)

cat("Sampling joint posterior from computed marginals...\n")
samples = inla.posterior.sample(2000, fit)

########################### SAVE RESULTS ###########################

# Only care about samples for betas (to save space)
beta_samples = t(vapply(samples, function(s) {
      s$latent[grep("^beta_", rownames(s$latent)), 1]
}, numeric(length(fit$summary.fixed$mean))))

colnames(beta_samples) = rownames(fit$summary.fixed)

# Full plot-ready marginals & joint posteriors for Beta; summary statistics for other marginals & hyperparameters
res = list(
  # store the controls used (so later scripts know what was requested)
  control = list(
    predictor = list(compute = TRUE),
    compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  ),

  # summaries
  beta_summary  = fit$summary.fixed,
  hyper_summary  = fit$summary.hyperpar,

  log_lambda = list(
      mean = fit$summary.linear.predictor$mean - log_size_factors,
      sd = fit$summary.linear.predictor$sd
  ),
  lambda = list(
      mean = fit$summary.fitted.values$mean / size_factors,
      sd = fit$summary.fitted.values$sd / size_factors
  ),
  epsilon = list(
      mean = fit$summary.random$sample_id$mean,
      sd = fit$summary.random$sample_id$sd
  ),

  # keep only beta marginals & samples to keep file size down
  beta_marginals = fit$marginals.fixed[colnames(X)],
  beta_samples = beta_samples
)
if (spatial) {
    res$psi = list(
        mean = fit$summary.random$region_id$mean,
        sd = fit$summary.random$region_id$sd
    )
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

cat("Saved results to:", output_file, "\n")
cat("Wrote file size (bytes):", file.info(output_file)$size, "\n")