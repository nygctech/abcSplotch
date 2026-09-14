#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript compute_contrasts.r <results_rds> <output_prefix> <delta_or_NA> <use_samples_TRUE_FALSE> [--beta-set=effective|components] <beta1:beta2> [beta3:beta4 ...]",
  "",
  "Arguments:",
  "  <results_rds>            Path to saved results object (.rds)",
  "  <output_prefix>          Prefix for output files",
  "  <delta_or_NA>            Minimum effect threshold (numeric) or NA",
  "  <use_samples_TRUE_FALSE> Whether to compute contrasts from stored posterior samples",
  "  --beta-set               'effective' (default) or 'components'",
  "  <beta1:beta2>            Contrast specification; names may come from any stored hierarchy level",
  "",
  "Examples:",
  "  Rscript compute_contrasts_hierarchical.r results.rds contrast_eval 2 FALSE beta_l2_c1_m6:beta_l2_c1_m2",
  "  Rscript compute_contrasts_hierarchical.r results.rds contrast_eval 2 TRUE --beta-set=effective beta_l3_c8_m6:beta_l3_c8_m2",
  "  Rscript compute_contrasts_hierarchical.r results.rds delta_eval NA TRUE --beta-set=components delta_l2_c8_m6:delta_l2_c8_m2",
  sep = "\n"
)

if (length(args) < 5L) stop(usage)

results_path <- args[1]
output_prefix <- args[2]
delta_arg <- args[3]
use_samples <- tolower(args[4]) %in% c("true", "t", "1", "yes")
remaining_args <- args[5:length(args)]

beta_set_arg <- grep("^--beta-set=", remaining_args, value = TRUE)
if (length(beta_set_arg) > 1L) stop("Specify --beta-set at most once")
beta_set <- if (length(beta_set_arg) == 1L) {
  sub("^--beta-set=", "", beta_set_arg)
} else {
  "effective"
}
remaining_args <- remaining_args[!grepl("^--beta-set=", remaining_args)]
contrast_specs <- remaining_args

if (!beta_set %in% c("effective", "components")) {
  stop("--beta-set must be 'effective' or 'components'")
}
if (length(contrast_specs) == 0L) stop("At least one contrast specification is required\n", usage)
if (!file.exists(results_path)) stop("Results file does not exist: ", results_path)

delta <- if (toupper(delta_arg) == "NA") NA_real_ else suppressWarnings(as.numeric(delta_arg))
if (!is.na(delta) && !is.finite(delta)) stop("delta must be numeric or NA")

res <- readRDS(results_path)

# Flatten a named list of per-level summary tables into one table.
flatten_beta_summaries <- function(summary_list) {
  if (is.null(summary_list) || !is.list(summary_list)) {
    stop("Expected a named list of beta summary tables")
  }

  out <- lapply(names(summary_list), function(level_name) {
    x <- summary_list[[level_name]]
    if (is.null(x)) return(NULL)
    if (!is.data.frame(x)) x <- as.data.frame(x)
    if (!"coef" %in% names(x)) {
      stop("Summary table '", level_name, "' does not contain a 'coef' column")
    }
    x$hierarchy_field <- level_name
    x
  })
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) stop("No beta summary tables were found")

  ret <- do.call(rbind, out)
  rownames(ret) <- NULL
  if (anyDuplicated(ret$coef)) {
    dup <- unique(ret$coef[duplicated(ret$coef)])
    stop("Coefficient names are duplicated across hierarchy fields: ", paste(dup, collapse = ", "))
  }
  ret
}

# Flatten a named list of sample matrices into one samples-by-coefficients matrix.
flatten_beta_samples <- function(sample_list) {
  if (is.null(sample_list) || !is.list(sample_list)) {
    stop("Expected a named list of beta sample matrices")
  }

  mats <- Filter(Negate(is.null), sample_list)
  if (length(mats) == 0L) stop("No beta sample matrices were found")

  n_rows <- vapply(mats, nrow, integer(1))
  if (length(unique(n_rows)) != 1L) {
    stop("All beta sample matrices must have the same number of posterior draws")
  }
  if (any(vapply(mats, function(x) is.null(colnames(x)), logical(1)))) {
    stop("All beta sample matrices must have coefficient column names")
  }

  ret <- do.call(cbind, mats)
  if (anyDuplicated(colnames(ret))) {
    dup <- unique(colnames(ret)[duplicated(colnames(ret))])
    stop("Coefficient names are duplicated across hierarchy fields: ", paste(dup, collapse = ", "))
  }
  ret
}

get_beta_summaries <- function(res, beta_set) {
  # New hierarchical structure.
  if (beta_set == "effective" && !is.null(res$beta_effective$summaries)) {
    return(flatten_beta_summaries(res$beta_effective$summaries))
  }
  if (beta_set == "components" && !is.null(res$beta_components$summaries)) {
    return(flatten_beta_summaries(res$beta_components$summaries))
  }

  # Legacy fallback applies only to effective coefficients.
  if (beta_set == "effective" && !is.null(res$beta_summary)) {
    warning("Using legacy res$beta_summary structure")
    return(res$beta_summary)
  }

  stop("Could not find ", beta_set, " beta summaries in the results object")
}

get_beta_samples <- function(res, beta_set) {
  # New hierarchical structure.
  if (beta_set == "effective" && !is.null(res$beta_effective$samples)) {
    return(flatten_beta_samples(res$beta_effective$samples))
  }
  if (beta_set == "components" && !is.null(res$beta_components$samples)) {
    return(flatten_beta_samples(res$beta_components$samples))
  }

  # Legacy fallback applies only to effective coefficients.
  if (beta_set == "effective" && !is.null(res$beta_samples)) {
    warning("Using legacy res$beta_samples structure")
    return(res$beta_samples)
  }

  stop(
    "Could not find ", beta_set, " beta samples in the results object. ",
    "The model must have been run with draw_samples=TRUE."
  )
}

compute_contrasts_from_marginals <- function(beta_summary, contrast_specs, delta = NA_real_) {
  out <- vector("list", length(contrast_specs))

  for (j in seq_along(contrast_specs)) {
    spec <- contrast_specs[j]
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) stop("Bad contrast spec: ", spec, ". Expected betaA:betaB")

    beta_plus <- trimws(parts[1])
    beta_minus <- trimws(parts[2])
    row_plus <- beta_summary[beta_summary$coef == beta_plus, , drop = FALSE]
    row_minus <- beta_summary[beta_summary$coef == beta_minus, , drop = FALSE]

    if (nrow(row_plus) == 0L) stop("Contrast beta not found: ", beta_plus)
    if (nrow(row_minus) == 0L) stop("Contrast beta not found: ", beta_minus)
    if (nrow(row_plus) > 1L) stop("Duplicated coefficient: ", beta_plus)
    if (nrow(row_minus) > 1L) stop("Duplicated coefficient: ", beta_minus)

    mu_delta <- row_plus$mean - row_minus$mean
    var_delta <- row_plus$sd^2 + row_minus$sd^2
    sd_delta <- sqrt(var_delta)
    p_gt0 <- 1 - pnorm(0, mean = mu_delta, sd = sd_delta)
    p_lt0 <- pnorm(0, mean = mu_delta, sd = sd_delta)

    row <- data.frame(
      contrast = paste0(beta_plus, "_minus_", beta_minus),
      plus = beta_plus,
      minus = beta_minus,
      beta_set = beta_set,
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
      row$p_gt_delta <- 1 - pnorm(delta, mean = mu_delta, sd = sd_delta)
      row$p_lt_neg_delta <- pnorm(-delta, mean = mu_delta, sd = sd_delta)
      row$p_large <- row$p_gt_delta + row$p_lt_neg_delta
      row$passes_prob_thresh <- pmax(p_gt0, p_lt0) > 0.95
      row$passes_effect_thresh <- abs(mu_delta) > delta
      row$passes_large_thresh <- row$p_large > 0.95
    }

    row$p_twosided <- 2 * min(p_gt0, p_lt0)
    row$p_tail <- if (!is.na(delta)) {
      pnorm(delta, mu_delta, sd_delta) - pnorm(-delta, mu_delta, sd_delta)
    } else NA_real_

    out[[j]] <- row
  }

  do.call(rbind, out)
}

compute_contrasts_from_samples <- function(beta_samples, contrast_specs, delta = NA_real_) {
  if (is.null(colnames(beta_samples))) stop("beta_samples must have column names")
  out <- vector("list", length(contrast_specs))

  for (j in seq_along(contrast_specs)) {
    spec <- contrast_specs[j]
    parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) stop("Bad contrast spec: ", spec, ". Expected betaA:betaB")

    beta_plus <- trimws(parts[1])
    beta_minus <- trimws(parts[2])
    if (!(beta_plus %in% colnames(beta_samples))) stop("Contrast beta not found: ", beta_plus)
    if (!(beta_minus %in% colnames(beta_samples))) stop("Contrast beta not found: ", beta_minus)

    d <- beta_samples[, beta_plus] - beta_samples[, beta_minus]
    mu_delta <- mean(d, na.rm = TRUE)
    var_delta <- stats::var(d, na.rm = TRUE)
    sd_delta <- sqrt(var_delta)
    p_gt0 <- mean(d > 0, na.rm = TRUE)
    p_lt0 <- mean(d < 0, na.rm = TRUE)

    row <- data.frame(
      contrast = paste0(beta_plus, "_minus_", beta_minus),
      plus = beta_plus,
      minus = beta_minus,
      beta_set = beta_set,
      source = "posterior_samples",
      mean = mu_delta,
      var = var_delta,
      sd = sd_delta,
      q025 = unname(stats::quantile(d, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(d, 0.975, na.rm = TRUE)),
      p_gt0 = p_gt0,
      p_lt0 = p_lt0,
      stringsAsFactors = FALSE
    )

    if (!is.na(delta)) {
      row$p_gt_delta <- mean(d > delta, na.rm = TRUE)
      row$p_lt_neg_delta <- mean(d < -delta, na.rm = TRUE)
      row$p_large <- mean(abs(d) > delta, na.rm = TRUE)
      row$passes_prob_thresh <- pmax(p_gt0, p_lt0) > 0.95
      row$passes_effect_thresh <- abs(mu_delta) > delta
      row$passes_large_thresh <- row$p_large > 0.95
    }

    row$p_twosided <- 2 * min(p_gt0, p_lt0)
    row$p_tail <- if (!is.na(delta)) mean(abs(d) <= delta, na.rm = TRUE) else NA_real_
    out[[j]] <- row
  }

  do.call(rbind, out)
}

if (use_samples) {
  beta_samples <- get_beta_samples(res, beta_set)
  contrast_tbl <- compute_contrasts_from_samples(beta_samples, contrast_specs, delta)
} else {
  beta_summary <- get_beta_summaries(res, beta_set)
  contrast_tbl <- compute_contrasts_from_marginals(beta_summary, contrast_specs, delta)
}

contrast_tbl$delta <- delta
saveRDS(contrast_tbl, paste0(output_prefix, ".rds"))
write.csv(contrast_tbl, paste0(output_prefix, ".csv"), row.names = FALSE)

cat("Saved contrast table to:\n")
cat("  ", paste0(output_prefix, ".rds"), "\n")
cat("  ", paste0(output_prefix, ".csv"), "\n")
print(contrast_tbl)
