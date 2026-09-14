#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript plot_beta_kde_hierarchical.R <results_rds> <mode> <beta_names_comma_sep> [labels_comma_sep] [output_png] [--beta-set=effective|components]",
  "",
  "Arguments:",
  "  <results_rds>          Path to saved RDS results object",
  "  <mode>                 Either 'marginal' or 'joint'",
  "  <beta_names_comma_sep> Comma-separated coefficient names",
  "  [labels_comma_sep]     Optional comma-separated legend labels; use NA to omit",
  "  [output_png]           Optional output PNG path; use NA for interactive plotting",
  "  --beta-set             'effective' (default) or 'components'",
  "",
  "Examples:",
  "  Rscript plot_beta_kde_hierarchical.R results.rds marginal beta_l2_c1_m1,beta_l2_c2_m1",
  "  Rscript plot_beta_kde_hierarchical.R results.rds joint beta_l3_c8_m1,beta_l3_c9_m1 Cond8,Cond9 plot.png",
  "  Rscript plot_beta_kde_hierarchical.R results.rds marginal delta_l2_c8_m1,delta_l2_c9_m1 NA delta.png --beta-set=components",
  sep = "\n"
)

if (length(args) < 3L) stop(usage)

beta_set_arg <- grep("^--beta-set=", args, value = TRUE)
if (length(beta_set_arg) > 1L) stop("Specify --beta-set at most once")
beta_set <- if (length(beta_set_arg) == 1L) sub("^--beta-set=", "", beta_set_arg) else "effective"
args <- args[!grepl("^--beta-set=", args)]

if (!beta_set %in% c("effective", "components")) {
  stop("--beta-set must be 'effective' or 'components'")
}
if (length(args) < 3L) stop(usage)

results_path <- args[1]
mode <- tolower(args[2])
beta_names <- trimws(strsplit(args[3], ",", fixed = TRUE)[[1]])
labels <- NULL
output_png <- NULL

if (length(args) >= 4L && toupper(args[4]) != "NA") {
  labels <- trimws(strsplit(args[4], ",", fixed = TRUE)[[1]])
}
if (length(args) >= 5L && toupper(args[5]) != "NA") output_png <- args[5]

if (!mode %in% c("marginal", "joint")) stop("mode must be either 'marginal' or 'joint'")
if (!file.exists(results_path)) stop("Results file does not exist: ", results_path)

res <- readRDS(results_path)
if (is.null(labels)) labels <- beta_names
if (length(labels) != length(beta_names)) stop("Number of labels must match number of beta names")

flatten_beta_marginals <- function(marginal_list) {
  if (is.null(marginal_list) || !is.list(marginal_list)) {
    stop("Expected a named list of beta marginal lists")
  }
  pieces <- lapply(marginal_list, function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
      stop("Every stored beta marginal must have a coefficient name")
    }
    x
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0L) stop("No beta marginals were found")

  # Concatenate without prefixing coefficient names by hierarchy-list names.
  out <- do.call(c, unname(pieces))
  if (anyDuplicated(names(out))) {
    dup <- unique(names(out)[duplicated(names(out))])
    stop("Coefficient names are duplicated across hierarchy fields: ", paste(dup, collapse = ", "))
  }
  out
}

flatten_beta_samples <- function(sample_list) {
  if (is.null(sample_list) || !is.list(sample_list)) {
    stop("Expected a named list of beta sample matrices")
  }
  mats <- Filter(Negate(is.null), sample_list)
  if (length(mats) == 0L) stop("No beta sample matrices were found")
  n_rows <- vapply(mats, nrow, integer(1))
  if (length(unique(n_rows)) != 1L) stop("All sample matrices must contain the same number of draws")
  if (any(vapply(mats, function(x) is.null(colnames(x)), logical(1)))) {
    stop("All sample matrices must have coefficient column names")
  }
  out <- do.call(cbind, mats)
  if (anyDuplicated(colnames(out))) {
    dup <- unique(colnames(out)[duplicated(colnames(out))])
    stop("Coefficient names are duplicated across hierarchy fields: ", paste(dup, collapse = ", "))
  }
  out
}

get_beta_marginals <- function(res, beta_set) {
  if (beta_set == "effective" && !is.null(res$beta_effective$marginals)) {
    return(flatten_beta_marginals(res$beta_effective$marginals))
  }
  if (beta_set == "components" && !is.null(res$beta_components$marginals)) {
    return(flatten_beta_marginals(res$beta_components$marginals))
  }
  if (beta_set == "effective" && !is.null(res$beta_marginals)) {
    warning("Using legacy res$beta_marginals structure")
    return(res$beta_marginals)
  }
  stop("Could not find ", beta_set, " beta marginals in the results object")
}

get_beta_samples <- function(res, beta_set) {
  if (beta_set == "effective" && !is.null(res$beta_effective$samples)) {
    return(flatten_beta_samples(res$beta_effective$samples))
  }
  if (beta_set == "components" && !is.null(res$beta_components$samples)) {
    return(flatten_beta_samples(res$beta_components$samples))
  }
  if (beta_set == "effective" && !is.null(res$beta_samples)) {
    warning("Using legacy res$beta_samples structure")
    return(res$beta_samples)
  }
  stop(
    "Could not find ", beta_set, " beta samples in the results object. ",
    "The model must have been run with draw_samples=TRUE."
  )
}

cols <- grDevices::hcl.colors(length(beta_names), palette = "Dark 3")
plot_title_suffix <- if (beta_set == "effective") "effective beta" else "hierarchy components"

if (mode == "marginal") {
  all_marginals <- get_beta_marginals(res, beta_set)
  missing <- setdiff(beta_names, names(all_marginals))
  if (length(missing) > 0L) stop("Missing beta marginals: ", paste(missing, collapse = ", "))
  marginals <- all_marginals[beta_names]

  valid_marginal <- vapply(
    marginals,
    function(m) is.matrix(m) && ncol(m) >= 2L && nrow(m) > 0L,
    logical(1)
  )
  if (!all(valid_marginal)) {
    stop("Invalid marginal matrices: ", paste(names(marginals)[!valid_marginal], collapse = ", "))
  }

  xlim <- range(unlist(lapply(marginals, function(m) m[, 1])), finite = TRUE)
  ylim <- c(0, max(unlist(lapply(marginals, function(m) m[, 2])), na.rm = TRUE))

  if (!is.null(output_png)) {
    grDevices::png(output_png, width = 900, height = 700, type = "cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  plot(
    marginals[[1]][, 1], marginals[[1]][, 2],
    type = "l", lwd = 2, col = cols[1], xlim = xlim, ylim = ylim,
    xlab = expression(beta), ylab = "Density",
    main = paste("Posterior marginals:", plot_title_suffix)
  )
  if (length(marginals) > 1L) {
    for (i in 2:length(marginals)) {
      lines(marginals[[i]][, 1], marginals[[i]][, 2], lwd = 2, col = cols[i])
    }
  }
  legend("topright", legend = labels, col = cols, lwd = 2, bty = "n")
}

if (mode == "joint") {
  beta_samples <- get_beta_samples(res, beta_set)
  missing <- setdiff(beta_names, colnames(beta_samples))
  if (length(missing) > 0L) stop("Missing beta samples: ", paste(missing, collapse = ", "))

  densities <- lapply(beta_names, function(b) {
    x <- beta_samples[, b]
    x <- x[is.finite(x)]
    if (length(x) < 2L) stop("Not enough finite posterior draws for ", b)
    stats::density(x)
  })
  names(densities) <- beta_names

  xlim <- range(unlist(lapply(densities, function(d) d$x)), finite = TRUE)
  ylim <- c(0, max(unlist(lapply(densities, function(d) d$y)), na.rm = TRUE))

  if (!is.null(output_png)) {
    grDevices::png(output_png, width = 900, height = 700, type = "cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  plot(
    densities[[1]]$x, densities[[1]]$y,
    type = "l", lwd = 2, col = cols[1], xlim = xlim, ylim = ylim,
    xlab = expression(beta), ylab = "Density",
    main = paste("KDE of joint posterior samples:", plot_title_suffix)
  )
  if (length(densities) > 1L) {
    for (i in 2:length(densities)) {
      lines(densities[[i]]$x, densities[[i]]$y, lwd = 2, col = cols[i])
    }
  }
  legend("topright", legend = labels, col = cols, lwd = 2, bty = "n")
}
