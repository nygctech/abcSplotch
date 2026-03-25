#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript plot_beta_kde.R <results_rds> <mode> <beta_names_comma_sep> [labels_comma_sep] [output_png]",
  "",
  "Arguments:",
  "  <results_rds>         Path to saved RDS results object",
  "  <mode>                Either 'marginal' or 'joint'",
  "  <beta_names_comma_sep> Comma-separated beta names, e.g. beta_c1_m1_k1,beta_c2_m1_k1",
  "  [labels_comma_sep]    Optional comma-separated legend labels",
  "  [output_png]          Optional output PNG path; if omitted, plot is shown interactively",
  "",
  "Examples:",
  "  Rscript plot_beta_kde.R results_1.R marginal beta_c1_m1_k1,beta_c2_m1_k1",
  "  Rscript plot_beta_kde.R results_1.R joint beta_c1_m1_k1,beta_c2_m1_k1 Cond1,Cond2 plot.png",
  sep = "\n"
)

if (length(args) < 3) {
  stop(usage)
}

results_path <- args[1]
mode <- tolower(args[2])
beta_names <- strsplit(args[3], ",", fixed = TRUE)[[1]]
beta_names <- trimws(beta_names)

labels <- NULL
output_png <- NULL

if (length(args) >= 4) {
  labels <- strsplit(args[4], ",", fixed = TRUE)[[1]]
  labels <- trimws(labels)
}

if (length(args) >= 5) {
  output_png <- args[5]
}

if (!mode %in% c("marginal", "joint")) {
  stop("mode must be either 'marginal' or 'joint'")
}

if (!file.exists(results_path)) {
  stop("Results file does not exist: ", results_path)
}

res <- readRDS(results_path)

if (is.null(labels)) {
  labels <- beta_names
}

if (length(labels) != length(beta_names)) {
  stop("Number of labels must match number of beta names")
}

cols <- grDevices::hcl.colors(length(beta_names), palette = "Dark 3")

# ---------------------------
# Marginal mode
# ---------------------------
if (mode == "marginal") {
  if (is.null(res$beta_marginals)) {
    stop("res$beta_marginals not found in RDS file")
  }

  missing <- setdiff(beta_names, names(res$beta_marginals))
  if (length(missing) > 0) {
    stop("Missing beta marginals: ", paste(missing, collapse = ", "))
  }

  marginals <- res$beta_marginals[beta_names]

  # Determine common plot limits
  xlim <- range(unlist(lapply(marginals, function(m) m[, 1])))
  ylim <- range(unlist(lapply(marginals, function(m) m[, 2])))

  if (!is.null(output_png)) {
    png(output_png, width = 900, height = 700, type = "cairo")
    on.exit(dev.off(), add = TRUE)
  }

  plot(
    marginals[[1]][, 1], marginals[[1]][, 2],
    type = "l",
    lwd = 2,
    col = cols[1],
    xlim = xlim,
    ylim = ylim,
    xlab = expression(beta),
    ylab = "Density",
    main = "Posterior marginals of beta"
  )

  if (length(marginals) > 1) {
    for (i in 2:length(marginals)) {
      lines(
        marginals[[i]][, 1],
        marginals[[i]][, 2],
        lwd = 2,
        col = cols[i]
      )
    }
  }

  legend("topright", legend = labels, col = cols, lwd = 2, bty = "n")
}

# ---------------------------
# Joint mode
# ---------------------------
if (mode == "joint") {
  if (is.null(res$beta_samples)) {
    stop("res$beta_samples not found in RDS file")
  }

  if (is.null(colnames(res$beta_samples))) {
    stop("res$beta_samples has no column names")
  }

  missing <- setdiff(beta_names, colnames(res$beta_samples))
  if (length(missing) > 0) {
    stop("Missing beta samples: ", paste(missing, collapse = ", "))
  }

  densities <- lapply(beta_names, function(b) density(res$beta_samples[, b], na.rm = TRUE))
  names(densities) <- beta_names

  xlim <- range(unlist(lapply(densities, function(d) d$x)))
  ylim <- range(unlist(lapply(densities, function(d) d$y)))

  if (!is.null(output_png)) {
    png(output_png, width = 900, height = 700, type = "cairo")
    on.exit(dev.off(), add = TRUE)
  }

  plot(
    densities[[1]]$x, densities[[1]]$y,
    type = "l",
    lwd = 2,
    col = cols[1],
    xlim = xlim,
    ylim = ylim,
    xlab = expression(beta),
    ylab = "Density",
    main = "KDE of joint posterior samples for beta"
  )

  if (length(densities) > 1) {
    for (i in 2:length(densities)) {
      lines(
        densities[[i]]$x,
        densities[[i]]$y,
        lwd = 2,
        col = cols[i]
      )
    }
  }

  legend("topright", legend = labels, col = cols, lwd = 2, bty = "n")
}