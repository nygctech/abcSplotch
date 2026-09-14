extract_latent_component <- function(
    posterior_sample,
    latent_name,
    n_coefficients
) {
    latent <- posterior_sample$latent

    if (is.matrix(latent)) {
        latent_values <- latent[, 1]
        latent_names <- rownames(latent)
    } else {
        latent_values <- as.numeric(latent)
        latent_names <- names(latent)
    }

    candidate_names <- list(
        paste0(latent_name, ":", seq_len(n_coefficients)),
        paste0(latent_name, seq_len(n_coefficients))
    )

    idx <- NULL

    for (candidate in candidate_names) {
        candidate_idx <- match(candidate, latent_names)

        if (!anyNA(candidate_idx)) {
            idx <- candidate_idx
            break
        }
    }

    # Fallback for names such as beta_l1:1 or beta_l1[1].
    if (is.null(idx)) {
        prefix_idx <- grep(
            paste0("^", latent_name, "(:|\\[|$)"),
            latent_names
        )

        if (length(prefix_idx) == n_coefficients) {
            idx <- prefix_idx
        }
    }

    if (is.null(idx) || length(idx) != n_coefficients) {
        stop(
            "Expected ",
            n_coefficients,
            " entries for latent effect '",
            latent_name,
            "', but could not identify them in the posterior sample."
        )
    }

    as.numeric(latent_values[idx])
}


sample_hierarchical_betas <- function(
    fit,
    beta_designs,
    effective_beta_metadata,
    n_hierarchy_levels,
    n_total = 1000L,
    chunk_size = 50L,
    beta_priors = FALSE,
    beta_mean = NULL,
    verbose = TRUE
) {
    stopifnot(
        n_hierarchy_levels >= 1L,
        length(beta_designs) == n_hierarchy_levels,
        n_total >= 1L,
        chunk_size >= 1L
    )

    beta_latent_names <- vapply(
        beta_designs,
        function(design) design$latent_name,
        character(1)
    )

    beta_dimensions <- vapply(
        beta_designs,
        function(design) as.integer(design$p),
        integer(1)
    )

    names(beta_dimensions) <- beta_latent_names

    expected_latent_names <- c(
        "beta_l1",
        if (n_hierarchy_levels >= 2L) "delta_l2",
        if (n_hierarchy_levels >= 3L) "delta_l3"
    )

    if (!all(expected_latent_names %in% beta_latent_names)) {
        stop(
            "Expected latent effects were not found in beta_designs. ",
            "Expected: ",
            paste(expected_latent_names, collapse = ", "),
            ". Found: ",
            paste(beta_latent_names, collapse = ", ")
        )
    }

    if (verbose) {
        cat(
            "Sampling joint posterior for hierarchical beta coefficients...\n"
        )
    }

    # One samples-by-coefficients matrix for each raw latent component.
    component_samples <- lapply(
        beta_dimensions,
        function(p_level) {
            matrix(
                NA_real_,
                nrow = n_total,
                ncol = p_level
            )
        }
    )

    for (l in seq_along(beta_designs)) {
        latent_name <- beta_designs[[l]]$latent_name
        coefficient_names <- beta_designs[[l]]$coefficient_names

        if (length(coefficient_names) != beta_designs[[l]]$p) {
            stop(
                "The number of coefficient names does not equal p for ",
                latent_name,
                "."
            )
        }

        colnames(component_samples[[latent_name]]) <-
            coefficient_names
    }

    posterior_selection <- lapply(
        beta_dimensions,
        seq_len
    )

    starts <- seq.int(
        from = 1L,
        to = n_total,
        by = chunk_size
    )

    for (s in starts) {
        m <- min(
            chunk_size,
            n_total - s + 1L
        )

        row_idx <- s:(s + m - 1L)

        samp_chunk <- INLA::inla.posterior.sample(
            n = m,
            result = fit,
            selection = posterior_selection,
            add.names = TRUE
        )

        for (latent_name in beta_latent_names) {
            p_level <- beta_dimensions[[latent_name]]

            component_samples[[latent_name]][row_idx, ] <- t(
                vapply(
                    samp_chunk,
                    extract_latent_component,
                    numeric(p_level),
                    latent_name = latent_name,
                    n_coefficients = p_level
                )
            )
        }

        rm(samp_chunk)
        gc(verbose = FALSE)

        if (verbose) {
            cat(
                "  sampled ",
                min(s + m - 1L, n_total),
                " / ",
                n_total,
                "\n",
                sep = ""
            )
        }
    }

    # Effective coefficients are constructed from the same joint draws,
    # preserving posterior covariance among hierarchy components.
    effective_samples <- list(
        beta_l1 = component_samples$beta_l1
    )

    if (n_hierarchy_levels >= 2L) {
        map_l2 <- effective_beta_metadata[[2]]

        required_l2_columns <- c(
            "parent_l1_col",
            "delta_l2_col",
            "coefficient"
        )

        if (!all(required_l2_columns %in% names(map_l2))) {
            stop(
                "Level-2 effective-beta metadata must contain: ",
                paste(required_l2_columns, collapse = ", ")
            )
        }

        effective_samples$beta_l2 <-
            component_samples$beta_l1[
                ,
                map_l2$parent_l1_col,
                drop = FALSE
            ] +
            component_samples$delta_l2[
                ,
                map_l2$delta_l2_col,
                drop = FALSE
            ]

        colnames(effective_samples$beta_l2) <-
            map_l2$coefficient
    }

    if (n_hierarchy_levels >= 3L) {
        map_l3 <- effective_beta_metadata[[3]]

        required_l3_columns <- c(
            "parent_l1_col",
            "parent_l2_col",
            "delta_l3_col",
            "coefficient"
        )

        if (!all(required_l3_columns %in% names(map_l3))) {
            stop(
                "Level-3 effective-beta metadata must contain: ",
                paste(required_l3_columns, collapse = ", ")
            )
        }

        effective_samples$beta_l3 <-
            component_samples$beta_l1[
                ,
                map_l3$parent_l1_col,
                drop = FALSE
            ] +
            component_samples$delta_l2[
                ,
                map_l3$parent_l2_col,
                drop = FALSE
            ] +
            component_samples$delta_l3[
                ,
                map_l3$delta_l3_col,
                drop = FALSE
            ]

        colnames(effective_samples$beta_l3) <-
            map_l3$coefficient
    }

    # Restore nonzero level-1 prior means. The delta terms remain centered;
    # each effective lower-level coefficient receives its level-1 offset.
    if (beta_priors) {
        if (is.null(beta_mean)) {
            stop("beta_mean must be supplied when beta_priors = TRUE.")
        }

        if (length(beta_mean) != ncol(component_samples$beta_l1)) {
            stop(
                "length(beta_mean) = ",
                length(beta_mean),
                ", but beta_l1 has ",
                ncol(component_samples$beta_l1),
                " coefficients."
            )
        }

        component_samples$beta_l1 <- sweep(
            component_samples$beta_l1,
            MARGIN = 2L,
            STATS = beta_mean,
            FUN = "+"
        )

        effective_samples$beta_l1 <-
            component_samples$beta_l1

        if (n_hierarchy_levels >= 2L) {
            shift_l2 <- beta_mean[
                effective_beta_metadata[[2]]$parent_l1_col
            ]

            effective_samples$beta_l2 <- sweep(
                effective_samples$beta_l2,
                MARGIN = 2L,
                STATS = shift_l2,
                FUN = "+"
            )
        }

        if (n_hierarchy_levels >= 3L) {
            shift_l3 <- beta_mean[
                effective_beta_metadata[[3]]$parent_l1_col
            ]

            effective_samples$beta_l3 <- sweep(
                effective_samples$beta_l3,
                MARGIN = 2L,
                STATS = shift_l3,
                FUN = "+"
            )
        }
    }

    list(
        components = component_samples,
        effective = effective_samples,
        n_samples = n_total
    )
}