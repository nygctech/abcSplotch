#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript compute_contrasts_from_res.R <results_rds> <output_prefix> <delta_or_NA> <use_samples_TRUE_FALSE> <beta1:beta2> [beta3:beta4 ...]",
  "",
  "Arguments:",
  "  <results_rds>            Path to saved results object (.rds)",
  "  <output_prefix>          Prefix for output files",
  "  <delta_or_NA>            Minimum effect threshold (numeric) or NA",
  "  <use_samples_TRUE_FALSE> Whether to compute contrasts from stored beta_samples",
  "  <beta1:beta2>            Contrast specification",
  "",
  "Examples:",
  "  Rscript compute_contrasts_from_res.R results_real.rds contrast_eval 2 FALSE beta_c1_m6:beta_c1_m2",
  "  Rscript compute_contrasts_from_res.R results_real.rds contrast_eval 2 TRUE  beta_c1_m6:beta_c1_m2",
  sep = "\n"
)

if (length(args) < 5) {
  stop(usage)
}

results_path <- args[1]
output_prefix <- args[2]
delta_arg <- args[3]
use_samples <- tolower(args[4]) %in% c("true", "t", "1", "yes")
contrast_specs <- args[5:length(args)]

if (!file.exists(results_path)) {
  stop("Results file does not exist: ", results_path)
}

delta <- if (toupper(delta_arg) == "NA") NA_real_ else as.numeric(delta_arg)
if (!is.na(delta) && !is.finite(delta)) {
  stop("delta must be numeric or NA")
}

res <- readRDS(results_path)

# ---------------------------
# Helper: contrasts from marginal summaries
# ---------------------------
compute_contrasts_from_marginals <- function(beta_summary, contrast_specs, delta = NA_real_) {
  out <- vector("list", length(contrast_specs))

  for (j in seq_along(contrast_specs)) {
    spec <- contrast_specs[j]
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]

    if (length(parts) != 2) {
      stop("Bad contrast spec: ", spec, ". Expected betaA:betaB")
    }

    beta_plus  <- trimws(parts[1])
    beta_minus <- trimws(parts[2])

    row_plus  <- beta_summary[beta_summary$coef == beta_plus, , drop = FALSE]
    row_minus <- beta_summary[beta_summary$coef == beta_minus, , drop = FALSE]

    if (nrow(row_plus) == 0) stop("Contrast beta not found in beta_summary: ", beta_plus)
    if (nrow(row_minus) == 0) stop("Contrast beta not found in beta_summary: ", beta_minus)
    if (nrow(row_plus) > 1) stop("Duplicated coef entry in beta_summary: ", beta_plus)
    if (nrow(row_minus) > 1) stop("Duplicated coef entry in beta_summary: ", beta_minus)

    mu_delta <- row_plus$mean - row_minus$mean
    var_delta <- row_plus$sd^2 + row_minus$sd^2
    sd_delta <- sqrt(var_delta)

    p_gt0 <- 1 - pnorm(0, mean = mu_delta, sd = sd_delta)
    p_lt0 <- pnorm(0, mean = mu_delta, sd = sd_delta)

    row <- data.frame(
      contrast = paste0(beta_plus, "_minus_", beta_minus),
      plus = beta_plus,
      minus = beta_minus,
      source = "marginal_approx",
      mean = mu_delta,
      var = var_delta,
      sd = sd_delta,
      q025 = qnorm(0.025, mean = mu_delta, sd = sd_delta),
      q975 = qnorm(0.975, mean = mu_delta, sd = sd_delta),
      p_gt0 = p_gt0,
      p_lt0 = p_lt0,
      stringsAsFactors = FALSE
    )

    if (!is.na(delta)) {
      p_gt_delta <- 1 - pnorm(delta, mean = mu_delta, sd = sd_delta)
      p_lt_neg_delta <- pnorm(-delta, mean = mu_delta, sd = sd_delta)
      p_large <- p_gt_delta + p_lt_neg_delta

      row$p_gt_delta <- p_gt_delta
      row$p_lt_neg_delta <- p_lt_neg_delta
      row$p_large <- p_large
      row$passes_prob_thresh <- pmax(p_gt0, p_lt0) > 0.95
      row$passes_effect_thresh <- abs(mu_delta) > delta
      row$passes_large_thresh <- p_large > 0.95
    }

    # p-value analogues:
    row$p_twosided <- 2 * min(p_gt0, p_lt0)
    row$p_tail <- if (!is.na(delta)) { 
        pnorm(delta, mu_delta, sd_delta) - pnorm(-delta, mu_delta, sd_delta) 
    } else NA_real_

    out[[j]] <- row
  }

  do.call(rbind, out)
}

# ---------------------------
# Helper: contrasts from posterior samples
# ---------------------------
compute_contrasts_from_samples <- function(beta_samples, contrast_specs, delta = NA_real_) {
  if (is.null(colnames(beta_samples))) {
    stop("beta_samples must have column names")
  }

  out <- vector("list", length(contrast_specs))

  for (j in seq_along(contrast_specs)) {
    spec <- contrast_specs[j]
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]

    if (length(parts) != 2) {
      stop("Bad contrast spec: ", spec, ". Expected betaA:betaB")
    }

    beta_plus  <- trimws(parts[1])
    beta_minus <- trimws(parts[2])

    if (!(beta_plus %in% colnames(beta_samples))) {
      stop("Contrast beta not found in beta_samples: ", beta_plus)
    }
    if (!(beta_minus %in% colnames(beta_samples))) {
      stop("Contrast beta not found in beta_samples: ", beta_minus)
    }

    d <- beta_samples[, beta_plus] - beta_samples[, beta_minus]

    mu_delta <- mean(d)
    var_delta <- stats::var(d)
    sd_delta <- sqrt(var_delta)

    p_gt0 <- mean(d > 0)
    p_lt0 <- mean(d < 0)

    row <- data.frame(
      contrast = paste0(beta_plus, "_minus_", beta_minus),
      plus = beta_plus,
      minus = beta_minus,
      source = "posterior_samples",
      mean = mu_delta,
      var = var_delta,
      sd = sd_delta,
      q025 = unname(stats::quantile(d, 0.025)),
      q975 = unname(stats::quantile(d, 0.975)),
      p_gt0 = p_gt0,
      p_lt0 = p_lt0,
      stringsAsFactors = FALSE
    )

    if (!is.na(delta)) {
      p_gt_delta <- mean(d > delta)
      p_lt_neg_delta <- mean(d < -delta)
      p_large <- mean(abs(d) > delta)

      row$p_gt_delta <- p_gt_delta
      row$p_lt_neg_delta <- p_lt_neg_delta
      row$p_large <- p_large
      row$passes_prob_thresh <- pmax(p_gt0, p_lt0) > 0.95
      row$passes_effect_thresh <- abs(mu_delta) > delta
      row$passes_large_thresh <- p_large > 0.95
    }

    # p-value analogues:
    row$p_twosided <- 2 * min(mean(d > 0), mean(d < 0))
    row$p_tail <- if (!is.na(delta)) { mean(abs(d) <= delta) } else NA_real_

    out[[j]] <- row
  }

  do.call(rbind, out)
}

# ---------------------------
# Dispatch
# ---------------------------
if (use_samples) {
  if (is.null(res$beta_samples)) {
    stop("use_samples=TRUE requested, but res$beta_samples is missing")
  }
  contrast_tbl <- compute_contrasts_from_samples(res$beta_samples, contrast_specs, delta)
} else {
  if (is.null(res$beta_summary)) {
    stop("use_samples=FALSE requested, but res$beta_summary is missing")
  }
  contrast_tbl <- compute_contrasts_from_marginals(res$beta_summary, contrast_specs, delta)
}
contrast_tbl$delta <- delta

saveRDS(contrast_tbl, paste0(output_prefix, ".rds"))
write.csv(contrast_tbl, paste0(output_prefix, ".csv"), row.names = FALSE)

cat("Saved contrast table to:\n")
cat("  ", paste0(output_prefix, ".rds"), "\n")
cat("  ", paste0(output_prefix, ".csv"), "\n")
print(contrast_tbl)