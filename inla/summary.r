# Collect raw summaries and marginals by hierarchy level 
collect_raw_beta_results <- function(fit, beta_designs) {
    summaries <- list()
    marginals <- list()
    metadata <- list()

    for (l in seq_along(beta_designs)) {
        design <- beta_designs[[l]]
        latent_name <- design$latent_name
        coefficient_names <- design$coefficient_names

        summary_l <- fit$summary.random[[latent_name]]
        marginals_l <- fit$marginals.random[[latent_name]]

        stopifnot(
            nrow(summary_l) == length(coefficient_names),
            length(marginals_l) == length(coefficient_names)
        )

        summary_l$coef <- coefficient_names
        names(marginals_l) <- coefficient_names

        summaries[[latent_name]] <- summary_l
        marginals[[latent_name]] <- marginals_l
        metadata[[latent_name]] <- design$coefficient_metadata
    }

    list(
        summaries = summaries,
        marginals = marginals,
        metadata = metadata
    )
}

# Collect effective summaries and marginals
# i.e, derive level 2 and 3 betas from relevant level 1 betas and deltas
collect_effective_beta_results <- function(
    fit,
    effective_beta_metadata,
    n_hierarchy_levels
) {
    summaries <- list()
    marginals <- list()
    metadata <- list()

    # Level 1 effective coefficients are simply beta_l1.
    beta_l1_summary <- fit$summary.random$beta_l1
    beta_l1_marginals <- fit$marginals.random$beta_l1

    summaries$beta_l1 <- beta_l1_summary
    marginals$beta_l1 <- beta_l1_marginals

    if (n_hierarchy_levels >= 2L) {
        lincomb_names_l2 <- effective_beta_metadata[[2]]$coefficient

        idx_l2 <- match(
            lincomb_names_l2,
            rownames(fit$summary.lincomb.derived)
        )

        if (anyNA(idx_l2)) {
            stop("Some level-2 effective beta summaries were not found.")
        }

        summaries$beta_l2 <- fit$summary.lincomb.derived[
            idx_l2,
            ,
            drop = FALSE
        ]

        summaries$beta_l2$coef <- lincomb_names_l2

        marginals$beta_l2 <- fit$marginals.lincomb.derived[
            lincomb_names_l2
        ]

        names(marginals$beta_l2) <- lincomb_names_l2
        metadata$beta_l2 <- effective_beta_metadata[[2]]
    }

    if (n_hierarchy_levels >= 3L) {
        lincomb_names_l3 <- effective_beta_metadata[[3]]$coefficient

        idx_l3 <- match(
            lincomb_names_l3,
            rownames(fit$summary.lincomb.derived)
        )

        if (anyNA(idx_l3)) {
            stop("Some level-3 effective beta summaries were not found.")
        }

        summaries$beta_l3 <- fit$summary.lincomb.derived[
            idx_l3,
            ,
            drop = FALSE
        ]

        summaries$beta_l3$coef <- lincomb_names_l3

        marginals$beta_l3 <- fit$marginals.lincomb.derived[
            lincomb_names_l3
        ]

        names(marginals$beta_l3) <- lincomb_names_l3
        metadata$beta_l3 <- effective_beta_metadata[[3]]
    }

    list(
        summaries = summaries,
        marginals = marginals,
        metadata = metadata
    )
}

# Compression to help with output file size
compress_marginal <- function(marginal, max_points = 200L) {
    if (is.null(marginal) || nrow(marginal) <= max_points) {
        return(marginal)
    }

    idx <- unique(
        round(
            seq(
                from = 1,
                to = nrow(marginal),
                length.out = max_points
            )
        )
    )

    marginal[idx, , drop = FALSE]
}

compress_marginal_list <- function(marginal_list, max_points = 200L) {
    lapply(
        marginal_list,
        compress_marginal,
        max_points = max_points
    )
}

# Applying shifts for nonzero beta prior means
shift_summary <- function(summary_df, shift) {
    stopifnot(nrow(summary_df) == length(shift))

    location_columns <- intersect(
        c(
            "mean",
            "mode",
            "0.025quant",
            "0.5quant",
            "0.975quant"
        ),
        names(summary_df)
    )

    for (column in location_columns) {
        summary_df[[column]] <- summary_df[[column]] + shift
    }

    summary_df
}

shift_marginals <- function(marginals, shift) {
    stopifnot(length(marginals) == length(shift))

    Map(
        function(marginal, value) {
            INLA::inla.tmarginal(
                fun = function(x) x + value,
                marginal = marginal
            )
        },
        marginals,
        shift
    )
}