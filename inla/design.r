# Construct "X" (condition x samples) design matrices for each level of
# the hierarchical sample annotation.
make_level_design <- function(
    level_labels,
    D,
    E,
    compositional,
    K,
    N,
    hierarchy_level
) {
    stopifnot(length(level_labels) == N)
    stopifnot(length(D) == N)

    level_fac <- factor(level_labels)
    mroi_fac  <- factor(D)

    # Use a separator unlikely to occur in the original labels.
    g_fac <- interaction(
        level_fac,
        mroi_fac,
        drop = TRUE,
        sep = ":"
    )

    g_id <- as.integer(g_fac)
    G <- nlevels(g_fac)

    if (compositional) {
        stopifnot(nrow(E) == N, ncol(E) == K)

        # For row i, assign E[i, 1:K] to the K columns belonging
        # to its hierarchy-level × MROI group.
        i_idx <- rep(seq_len(N), each = K)
        j_idx <- (
            rep(g_id, each = K) - 1L
        ) * K + rep(seq_len(K), times = N)

        x_val <- as.numeric(t(E))

        X <- Matrix::sparseMatrix(
            i = i_idx,
            j = j_idx,
            x = x_val,
            dims = c(N, G * K),
            giveCsparse = TRUE
        )

    } else {
        X <- Matrix::sparseMatrix(
            i = seq_len(N),
            j = g_id,
            x = 1,
            dims = c(N, G),
            giveCsparse = TRUE
        )
    }

    X <- Matrix::drop0(X)

    # Parse each observed hierarchy-level × MROI combination.
    g_levels <- levels(g_fac)
    parts <- strsplit(g_levels, ":", fixed = TRUE)

    level_label <- vapply(parts, `[`, character(1), 1)
    mroi_label  <- vapply(parts, `[`, character(1), 2)

    stopifnot(
        length(level_label) == G,
        length(mroi_label) == G
    )

    # Level 1 contains the parent beta terms.
    # Levels 2+ contain deviations from the preceding hierarchy.
    parameter_prefix <- if (hierarchy_level == 1L) {
        "beta_l1"
    } else {
        paste0("delta_l", hierarchy_level)
    }

    if (compositional) {
        coef_names <- unlist(
            lapply(seq_len(G), function(g) {
                paste0(
                    parameter_prefix,
                    "_c", level_label[g],
                    "_m", mroi_label[g],
                    "_k", seq_len(K)
                )
            }),
            use.names = FALSE
        )
    } else {
        coef_names <- vapply(
            seq_len(G),
            function(g) {
                paste0(
                    parameter_prefix,
                    "_c", level_label[g],
                    "_m", mroi_label[g]
                )
            },
            character(1)
        )
    }

    stopifnot(length(coef_names) == ncol(X))
    colnames(X) <- coef_names

    # In a compositional model, columns may be all zero if a cell type
    # has zero weight for every observation in a given group.
    nonzero_cols <- Matrix::colSums(abs(X)) > 0

    cat(
        "Hierarchy level", hierarchy_level, ": dropping",
        sum(!nonzero_cols), "all-zero coefficient columns\n"
    )

    X <- X[, nonzero_cols, drop = FALSE]
    X <- Matrix::drop0(X)

    # Retain metadata describing every fitted coefficient.
    coef_metadata <- if (compositional) {
        do.call(
            rbind,
            lapply(seq_len(G), function(g) {
                data.frame(
                    hierarchy_level = hierarchy_level,
                    condition = level_label[g],
                    mroi = mroi_label[g],
                    cell_type = seq_len(K),
                    coefficient = paste0(
                        parameter_prefix,
                        "_c", level_label[g],
                        "_m", mroi_label[g],
                        "_k", seq_len(K)
                    ),
                    stringsAsFactors = FALSE
                )
            })
        )
    } else {
        data.frame(
            hierarchy_level = hierarchy_level,
            condition = level_label,
            mroi = mroi_label,
            coefficient = coef_names,
            stringsAsFactors = FALSE
        )
    }

    coef_metadata <- coef_metadata[nonzero_cols, , drop = FALSE]
    rownames(coef_metadata) <- NULL

    list(
        X = X,
        p = ncol(X),
        hierarchy_level = hierarchy_level,
        latent_name = if (hierarchy_level == 1L) {
            "beta_l1"
        } else {
            paste0("delta_l", hierarchy_level)
        },
        group_factor = g_fac,
        group_levels = g_levels,
        level_labels = level_label,
        mroi_labels = mroi_label,
        retained_columns = nonzero_cols,
        coefficient_names = colnames(X),
        coefficient_metadata = coef_metadata
    )
}

# Construct a matching key for a coefficient while preserving MROI
# and compositional cell-type identity.
make_coef_key <- function(condition, mroi, cell_type = NULL) {
    if (is.null(cell_type)) {
        paste(condition, mroi, sep = "|")
    } else {
        paste(condition, mroi, cell_type, sep = "|")
    }
}

get_metadata_keys <- function(metadata) {
    if ("cell_type" %in% names(metadata)) {
        make_coef_key(
            condition = metadata$condition,
            mroi = metadata$mroi,
            cell_type = metadata$cell_type
        )
    } else {
        make_coef_key(
            condition = metadata$condition,
            mroi = metadata$mroi
        )
    }
}

# Level 2 and 3 shared effects are modeled as deviations from a level 1 term
# i.e., beta_l2 = beta_l1 + delta_2; beta_l3 = beta_l1 + delta_2 + delta_3
# This reconstructs the "effective" level 2 and 3 betas for saving & DE testing
make_effective_beta_lincombs <- function(
    beta_designs,
    rdat,
    n_hierarchy_levels
) {
    stopifnot(n_hierarchy_levels >= 1L)

    meta_l1 <- beta_designs[[1]]$coefficient_metadata
    keys_l1 <- get_metadata_keys(meta_l1)

    p_l1 <- beta_designs[[1]]$p

    lincombs <- list()
    lincomb_metadata <- list()

    # ----------------------------------------------------------
    # Level 1: no lincomb is required because beta_l1 is already
    # directly represented in the latent field.
    # ----------------------------------------------------------

    # ----------------------------------------------------------
    # Level 2 effective coefficients:
    # beta_l2 = beta_l1[parent] + delta_l2
    # ----------------------------------------------------------
    if (n_hierarchy_levels >= 2L) {
        meta_l2 <- beta_designs[[2]]$coefficient_metadata
        p_l2 <- beta_designs[[2]]$p
    
        parent_l1_condition <- rdat$level_2_mapping[
            as.integer(meta_l2$condition)
        ]
    
        if ("cell_type" %in% names(meta_l2)) {
            parent_keys <- make_coef_key(
                condition = parent_l1_condition,
                mroi = meta_l2$mroi,
                cell_type = meta_l2$cell_type
            )
        } else {
            parent_keys <- make_coef_key(
                condition = parent_l1_condition,
                mroi = meta_l2$mroi
            )
        }
    
        parent_l1_col <- match(parent_keys, keys_l1)
    
        if (anyNA(parent_l1_col)) {
            stop(
                "Could not match ",
                sum(is.na(parent_l1_col)),
                " level-2 coefficients to level-1 parents."
            )
        }
    
        for (j in seq_len(p_l2)) {
            w_l1 <- numeric(p_l1)
            w_l2 <- numeric(p_l2)
    
            w_l1[parent_l1_col[j]] <- 1
            w_l2[j] <- 1
    
            effective_name <- sub(
                "^delta_l2",
                "beta_l2",
                meta_l2$coefficient[j]
            )
    
            lc <- INLA::inla.make.lincomb(
                beta_l1 = w_l1,
                delta_l2 = w_l2
            )
    
            names(lc) <- effective_name
            lincombs <- c(lincombs, lc)
        }
    
        lincomb_metadata[[2]] <- data.frame(
            hierarchy_level = 2L,
            coefficient = sub(
                "^delta_l2",
                "beta_l2",
                meta_l2$coefficient
            ),
            level_1_condition = parent_l1_condition,
            level_2_condition = meta_l2$condition,
            mroi = meta_l2$mroi,
            parent_l1_col = parent_l1_col,
            delta_l2_col = seq_len(p_l2),
            stringsAsFactors = FALSE
        )
    
        if ("cell_type" %in% names(meta_l2)) {
            lincomb_metadata[[2]]$cell_type <- meta_l2$cell_type
        }
    }

    # ----------------------------------------------------------
    # Level 3 effective coefficients:
    # beta_l3 = beta_l1[parent] +
    #           delta_l2[parent] +
    #           delta_l3
    # ----------------------------------------------------------
    if (n_hierarchy_levels >= 3L) {
        meta_l2 <- beta_designs[[2]]$coefficient_metadata
        meta_l3 <- beta_designs[[3]]$coefficient_metadata
    
        keys_l2 <- get_metadata_keys(meta_l2)
    
        p_l2 <- beta_designs[[2]]$p
        p_l3 <- beta_designs[[3]]$p
    
        parent_l2_condition <- rdat$level_3_mapping[
            as.integer(meta_l3$condition)
        ]
    
        parent_l1_condition <- rdat$level_2_mapping[
            as.integer(parent_l2_condition)
        ]
    
        if ("cell_type" %in% names(meta_l3)) {
            parent_l2_keys <- make_coef_key(
                condition = parent_l2_condition,
                mroi = meta_l3$mroi,
                cell_type = meta_l3$cell_type
            )
    
            parent_l1_keys <- make_coef_key(
                condition = parent_l1_condition,
                mroi = meta_l3$mroi,
                cell_type = meta_l3$cell_type
            )
        } else {
            parent_l2_keys <- make_coef_key(
                condition = parent_l2_condition,
                mroi = meta_l3$mroi
            )
    
            parent_l1_keys <- make_coef_key(
                condition = parent_l1_condition,
                mroi = meta_l3$mroi
            )
        }
    
        parent_l2_col <- match(parent_l2_keys, keys_l2)
        parent_l1_col <- match(parent_l1_keys, keys_l1)
    
        if (anyNA(parent_l1_col)) {
            stop(
                "Could not match ",
                sum(is.na(parent_l1_col)),
                " level-3 coefficients to level-1 ancestors."
            )
        }
    
        if (anyNA(parent_l2_col)) {
            stop(
                "Could not match ",
                sum(is.na(parent_l2_col)),
                " level-3 coefficients to level-2 parents."
            )
        }
    
        for (j in seq_len(p_l3)) {
            w_l1 <- numeric(p_l1)
            w_l2 <- numeric(p_l2)
            w_l3 <- numeric(p_l3)
    
            w_l1[parent_l1_col[j]] <- 1
            w_l2[parent_l2_col[j]] <- 1
            w_l3[j] <- 1
    
            effective_name <- sub(
                "^delta_l3",
                "beta_l3",
                meta_l3$coefficient[j]
            )
    
            lc <- INLA::inla.make.lincomb(
                beta_l1 = w_l1,
                delta_l2 = w_l2,
                delta_l3 = w_l3
            )
    
            names(lc) <- effective_name
            lincombs <- c(lincombs, lc)
        }
    
        lincomb_metadata[[3]] <- data.frame(
            hierarchy_level = 3L,
            coefficient = sub(
                "^delta_l3",
                "beta_l3",
                meta_l3$coefficient
            ),
            level_1_condition = parent_l1_condition,
            level_2_condition = parent_l2_condition,
            level_3_condition = meta_l3$condition,
            mroi = meta_l3$mroi,
            parent_l1_col = parent_l1_col,
            parent_l2_col = parent_l2_col,
            delta_l3_col = seq_len(p_l3),
            stringsAsFactors = FALSE
        )
    
        if ("cell_type" %in% names(meta_l3)) {
            lincomb_metadata[[3]]$cell_type <- meta_l3$cell_type
        }
    }

    list(
        lincombs = lincombs,
        metadata = lincomb_metadata
    )
}