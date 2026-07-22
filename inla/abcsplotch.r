library(INLA)
library(rstan)
library(Matrix)
library(argparse)
source("design.r")
source("summary.r")
source("sampling.r")

script_start_time <- Sys.time()

########################### READ INPUT DATA ###########################

library(argparse)

parser <- ArgumentParser(
  description = "Run abcSplotch model with INLA inference."
)

# Positional arguments
parser$add_argument(
  "input_rdat_path",
  type = "character",
  help = "Path to input .rdat file"
)

parser$add_argument(
  "output_directory",
  type = "character",
  help = "Directory for model outputs"
)

# Optional arguments
parser$add_argument(
  "--draw-samples",
  action = "store_true",
  default = FALSE,
  help = "Draw posterior samples"
)

parser$add_argument(
  "--regional-precision",
  type = "double",
  default = log(200),
  help = "If provided, fixes BYM2 regional precision to this value (default: log(200))"
)

# Parse arguments
args <- parser$parse_args()

# Assign variables
input_path <- args$input_rdat_path
output_root <- args$output_directory
draw_samples <- args$draw_samples
regional_precision <- args$regional_precision

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
beta_priors = all(c("beta_prior_mean", "beta_prior_std") %in% names(rdat))
levels = rdat$N_levels

# Check validity of input
stopifnot(length(log_size_factors) == N)
stopifnot(length(D) == N)
stopifnot(length(y) == N)

# Check whether compositional and/or spatial information is present; adjust model accordingly
compositional = ("N_celltypes" %in% names(rdat))
spatial = (rdat$car > 0)
regional = ("region_list" %in% names(rdat))

if (compositional) {
    K = rdat$N_celltypes
    E = rdat$E
}

# Create sparse adjacency matrix
if (spatial) {
    # Map cells to modeled spatial regions, if desired
    if (regional) {
        region_id = rdat$region_list
    	R = max(region_id)

        stopifnot(length(region_id) == N)

        region_mapping = sparseMatrix(
            i = seq_len(N),
            j = region_id,
            x = 1,
            dims = c(N, R),
            giveCsparse = TRUE
        )
        W_dims = c(R, R)
    } else {
        W_dims = c(N, N)
    }
    i <- rdat$W_sparse[, 1]
    j <- rdat$W_sparse[, 2]
    
    # Build symmetric adjacency directly (still sparse CSC / dgCMatrix)
    W_sparse <- sparseMatrix(
      i = c(i, j),
      j = c(j, i),
      x = 1,
      dims = W_dims,
      giveCsparse = TRUE
    )
    
    # If duplicates exist, sparseMatrix will sum them (2,3,...) -> clamp back to 1
    W_sparse@x[] <- 1
    
    # Ensure no self edges (optional safety)
    Matrix::diag(W_sparse) <- 0
    W_sparse <- Matrix::drop0(W_sparse)
}

# Derive spot-level labeling for each covariate group
if ("N_level_3" %in% names(rdat) && rdat$N_level_3 > 0) {
    level_3 = rep(rdat$tissue_mapping, times = rdat$N_spots)
    level_2 = rdat$level_3_mapping[level_3]
    level_1 = rdat$level_2_mapping[level_2]
} else if ("N_level_2" %in% names(rdat) && rdat$N_level_2 > 0) {
    level_2 = rep(rdat$tissue_mapping, times = rdat$N_spots)
    level_1 = rdat$level_2_mapping[level_2]
} else {
    stopifnot("N_level_1" %in% names(rdat))
    stopifnot(rdat$N_level_1 > 0)
    level_1 = rep(rdat$tissue_mapping, times = rdat$N_spots)
}

########################### CONSTRUCT DESIGN MATRIX ###########################

beta_designs <- list()

beta_designs[[1]] <- make_level_design(
    level_labels = level_1,
    D = D,
    E = E,
    compositional = compositional,
    K = K,
    N = N,
    hierarchy_level = 1L
)

if (levels >= 2) {
    beta_designs[[2]] <- make_level_design(
        level_labels = level_2,
        D = D,
        E = E,
        compositional = compositional,
        K = K,
        N = N,
        hierarchy_level = 2L
    )
}

if (levels == 3) {
    beta_designs[[3]] <- make_level_design(
        level_labels = level_3,
        D = D,
        E = E,
        compositional = compositional,
        K = K,
        N = N,
        hierarchy_level = 3L
    )
}

n_hierarchy_levels <- length(beta_designs)

stopifnot(n_hierarchy_levels >= 1L)
stopifnot(
    all(vapply(
        beta_designs,
        function(x) nrow(x$X) == N,
        logical(1)
    ))
)

# Construct hierarchy-level beta/deviation effects and mappings.
effects <- list()
A <- list()

for (l in seq_along(beta_designs)) {
    design <- beta_designs[[l]]

    effects[[design$latent_name]] <- seq_len(design$p)
    A[[design$latent_name]] <- design$X
}

# Observation-level iid component.
effects$sample_id <- seq_len(N)
A$sample_id <- 1

# Optional spatial component.
if (spatial) {
    if (regional) {
        effects$region_id <- seq_len(R)
        A$region_id <- region_mapping
    } else {
        effects$region_id <- seq_len(N)
        A$region_id <- 1
    }
}

stopifnot(length(A) == length(effects))
stopifnot(identical(names(A), names(effects)))

stk <- inla.stack(
    data = list(y = y),
    A = A,
    effects = effects,
    tag = "est"
)

idx_est <- inla.stack.index(stk, tag = "est")$data

stopifnot(length(idx_est) == N)
stopifnot(length(idx_est) == length(size_factors))

cat("N:", N, "\n")
cat("levels:", levels, "\n")
cat("compositional:", compositional, "\n")
cat("spatial:", spatial, "\n")
if (spatial) {
    cat("regional:", regional, "\n")
    if (regional) { cat("R:", R, "\n") }
}
cat("length offset:", length(log_size_factors), "\n")
cat("nrow A:", nrow(inla.stack.A(stk)), "\n")
cat("length idx_est:", length(idx_est), "\n")

########################### DEFINE PRIORS ###########################

# BETA:
if (beta_priors) {
    mu = rdat$beta_prior_mean
    sigma = rdat$beta_prior_std
    
    if (compositional) {
        stopifnot(length(mu)==K)
        stopifnot(length(sigma)==K)
        b_idx = as.integer(sub(".*_k([0-9]+)$", "\\1", colnames(X)))
    } else {
        stopifnot(length(mu)==max(D))
        stopifnot(length(sigma)==max(D))
        b_idx = as.integer(sub(".*_m([0-9]+)$", "\\1", colnames(X)))
    }
    stopifnot(all(sigma > 0))
    stopifnot(!any(is.na(b_idx)))
    
    beta_mean = mu[b_idx]
    beta_prec = 1 / sigma[b_idx]^2

    Q_beta <- Matrix::Diagonal(
        n = p,
        x = beta_prec
    )
    beta_prior_offset <- as.vector(X %*% beta_mean)
} else {
    beta_prec = 0.25   # SD = 2, as in Stan formulation
}

# BETA LEVELS 2-3: 
# Additional variance on top of level 1 follows Half-Normal(0,1) in Splotch;
# this adapts for the log-precision speficiation employed by INLA
half_normal_sd1_prior <- "
expression:
    logdens = -0.5 * exp(-theta) - 0.5 * theta;
    return(logdens);
"

# EPSILON:
# In Stan implementation, hyperprior on SD of epsilon is HalfNormal(0,0.3)
# These PC prior settings make sure that SD of epsilon rarely exceeds what we'd expect from this
hyper_iid <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.773, 0.01)  # P(sigma > 0.773) = 0.01
  )
)

# PSI:
# LCAR only - for cell/spot-level cases where iid error added separately (epsilon)
hyper_besag = list(
    prec   = list(prior = "loggamma", param = c(1, 5e-4)),
    lambda = list(prior = "gaussian", param = c(0, 0.45))
)

# LCAR+iid - for region-level cases
hyper_bym2 <- list(
  # precision
  theta1 = list(
    initial = regional_precision,
    fixed = TRUE
  ),
  # phi (0=intraregion iid, 1=fully spatial)
  theta2 = list(
    prior = "pc",
    param = c(0.5, 0.8)
  )
)

########################### BUILD & RUN MODEL ###########################

# Formula: no intercept; offset passed later to fit; beta is a latent coefficient vector
if (beta_priors) {
    fml = y ~ 0 +
      f(beta_l1,
        model = "generic0",
        Cmatrix = Q_beta,
        constr = FALSE,
        hyper = list(
          prec = list(initial = 0, fixed = TRUE)
        )
      )    
} else {
    fml = y ~ 0 +
      f(beta_l1, model = "iid",
        hyper = list(
          prec = list(initial = log(beta_prec), fixed = TRUE)
        )
      )
}

# Additional hierarchical levels are supported through added variance ("delta") terms
if (levels >= 2L) {
    fml <- update(
        fml,
        . ~ . +
            f(
                delta_l2,
                model = "iid",
                hyper = half_normal_sd1_prior
            )
    )
}

if (levels >= 3L) {
    fml <- update(
        fml,
        . ~ . +
            f(
                delta_l3,
                model = "iid",
                hyper = half_normal_sd1_prior
            )
    )
}

# Define CAR component in one of two ways:
# 1. Regional: bym2 model mixes structured (inter-region) and unstructured (intra-region) with parameter "phi"
# 2. Node-level: besagproper2 pure Leroux-style structured variation dependent on neighbors
if (spatial) {
    if (regional) {
        fml = update(fml, . ~ . + f(region_id,
                                    model='bym2',
                                    graph=W_sparse,
                                    hyper=hyper_bym2))        
    } else {
        fml = update(fml, . ~ . + f(region_id,
                                    model='besagproper2',
                                    graph=W_sparse,
                                    hyper=hyper_besag))
    }
}
# For node-level case, iid epsilon term handles node/cell-level variation; in regional case this is overparameterized
# as the bym2 model allows for region-level iid effect
if (!regional) {
    fml = update(fml, . ~ . + + f(sample_id, model = "iid", hyper = hyper_iid))
}

# Compute "effective" level 2 and 3 betas by combining level 1 with modeled offsets
# (Need to specify that these combinations are modeled *before* fitting!)
effective_beta_def <- make_effective_beta_lincombs(
    beta_designs = beta_designs,
    rdat = rdat,
    n_hierarchy_levels = n_hierarchy_levels
)
effective_beta_lincombs <- effective_beta_def$lincombs
effective_beta_metadata <- effective_beta_def$metadata

fit_args = list(
    formula = fml,
    family = "nbinomial",
    data = inla.stack.data(stk),
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    offset = log_size_factors,
    lincomb = effective_beta_lincombs
)
if (beta_priors) { 
    fit_args$offset = log_size_factors + beta_prior_offset
}
fit = do.call(INLA::inla, fit_args)

# Quick check of level-1 beta terms
beta_l1_summary <- fit$summary.random$beta_l1
beta_l1_summary$coef <- colnames(beta_designs[[1]]$X)

stopifnot(nrow(beta_l1_summary) == ncol(beta_designs[[1]]$X))

print(beta_l1_summary[seq_len(min(10, nrow(beta_l1_summary))), ])
print(fit$summary.hyperpar)

########################### SAMPLE JOINT POSTERIOR ###########################

# Expensive; primarily retained for analyses requiring the covariance
# among hierarchical beta components
beta_posterior_samples <- NULL

if (draw_samples) {
    beta_posterior_samples <- sample_hierarchical_betas(
        fit = fit,
        beta_designs = beta_designs,
        effective_beta_metadata = effective_beta_metadata,
        n_hierarchy_levels = n_hierarchy_levels,
        n_total = 1000L,
        chunk_size = 50L,
        beta_priors = beta_priors,
        beta_mean = if (beta_priors) beta_mean else NULL
    )
}

########################### SAVE RESULTS ###########################

raw_beta <- collect_raw_beta_results(
    fit = fit,
    beta_designs = beta_designs
)

effective_beta <- collect_effective_beta_results(
    fit = fit,
    effective_beta_metadata = effective_beta_metadata,
    n_hierarchy_levels = n_hierarchy_levels
)

# Apply level-1 prior-mean offsets.
if (beta_priors) {
    raw_beta$summaries$beta_l1 <- shift_summary(
        raw_beta$summaries$beta_l1,
        beta_mean
    )

    raw_beta$marginals$beta_l1 <- shift_marginals(
        raw_beta$marginals$beta_l1,
        beta_mean
    )

    if (n_hierarchy_levels >= 2L) {
        shift_l2 <- beta_mean[
            effective_beta_metadata[[2]]$parent_l1_col
        ]

        effective_beta$summaries$beta_l2 <- shift_summary(
            effective_beta$summaries$beta_l2,
            shift_l2
        )

        effective_beta$marginals$beta_l2 <- shift_marginals(
            effective_beta$marginals$beta_l2,
            shift_l2
        )
    }

    if (n_hierarchy_levels >= 3L) {
        shift_l3 <- beta_mean[
            effective_beta_metadata[[3]]$parent_l1_col
        ]

        effective_beta$summaries$beta_l3 <- shift_summary(
            effective_beta$summaries$beta_l3,
            shift_l3
        )

        effective_beta$marginals$beta_l3 <- shift_marginals(
            effective_beta$marginals$beta_l3,
            shift_l3
        )
    }
}

# Level 1 effective beta is the adjusted raw level-1 beta.
effective_beta$summaries$beta_l1 <-
    raw_beta$summaries$beta_l1

effective_beta$marginals$beta_l1 <-
    raw_beta$marginals$beta_l1

effective_beta$metadata$beta_l1 <-
    raw_beta$metadata$beta_l1

# Reduce marginal storage size.
raw_beta$marginals <- lapply(
    raw_beta$marginals,
    compress_marginal_list,
    max_points = 200L
)

effective_beta$marginals <- lapply(
    effective_beta$marginals,
    compress_marginal_list,
    max_points = 200L
)

# ---- Returned structure ----
# -- Raw betas (for l2; deltas for l2-3) --
# res$beta_components$summaries
# res$beta_components$marginals
# res$beta_components$samples   (if joint samples drawn)
#
# -- Effective betas (l1-3) --
# res$beta_effective$summaries
# res$beta_effective$marginals
# res$beta_effective$samples   (if joint samples drawn)
res <- list(
    control = list(
        predictor = list(
            compute = TRUE
        ),
        compute = list(
            dic = TRUE,
            waic = TRUE,
            cpo = TRUE,
            config = TRUE
        )
    ),

    hierarchy = list(
        n_levels = n_hierarchy_levels,
        level_2_mapping = if (n_hierarchy_levels >= 2L) {
            rdat$level_2_mapping
        } else {
            NULL
        },
        level_3_mapping = if (n_hierarchy_levels >= 3L) {
            rdat$level_3_mapping
        } else {
            NULL
        }
    ),

    # Raw latent hierarchy:
    # beta_l1, delta_l2, delta_l3
    beta_components = list(
        summaries = raw_beta$summaries,
        marginals = raw_beta$marginals,
        metadata = raw_beta$metadata
    ),

    # Biologically interpreted coefficients:
    # beta_l1
    # beta_l1 + delta_l2
    # beta_l1 + delta_l2 + delta_l3
    beta_effective = list(
        summaries = effective_beta$summaries,
        marginals = effective_beta$marginals,
        metadata = effective_beta$metadata
    ),

    hyper_summary = fit$summary.hyperpar,

    log_lambda = list(
        mean = fit$summary.linear.predictor$mean[idx_est],
        sd = fit$summary.linear.predictor$sd[idx_est]
    ),

    lambda = list(
        mean = fit$summary.fitted.values$mean[idx_est] /
            size_factors,
        sd = fit$summary.fitted.values$sd[idx_est] /
            size_factors
    ),

    epsilon = list(
        mean = fit$summary.random$sample_id$mean,
        sd = fit$summary.random$sample_id$sd
    )
)

if (spatial) {
    res$psi <- list(
        mean = fit$summary.random$region_id$mean,
        sd = fit$summary.random$region_id$sd
    )
}

if (draw_samples) {
    res$beta_components$samples <-
        beta_posterior_samples$components

    res$beta_effective$samples <-
        beta_posterior_samples$effective
}      

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