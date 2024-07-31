#' Apply dimensionality reduction using WGCNA
#'
#' Performs Weighted gene correlation network analysis (WGCNA).
#' Calls \link{ReducedExperiment}[run_wgcna] to perform the analysis.
#'
#' @param X Either a \link{SummarizedExperiment}[SummarizedExperiment] object
#' or a matrix containing data to be subject to WGCNA. `X` should have rows as
#' features and columns as samples.
#'
#' @param center_X If TRUE, X is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to WGCNA.
#'
#' @param scale_X If TRUE, X is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before WGCNA.
#'
#' @param assay_name If `X` is a
#' \link{SummarizedExperiment}[SummarizedExperiment], then this should be the
#' name of the assay to be subject to WGCNA.
#'
#' @param ... Additional arguments to be passed to
#' \link{ReducedExperiment}[run_wgcna].
#'
#' @export
identify_modules <- function(
        X,
        center_X = TRUE,
        scale_X = TRUE,
        assay_name = "normal",
        ...) {
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    if (assay_name != "normal") {
        assay(X, "normal") <- assay(X, assay_name)
    }

    if ("transformed" %in% assayNames(X)) {
        warning("Overwriting 'transformed' assay slot in X")
    }
    assay(X, "transformed") <-
        t(scale(t(assay(X, "normal")), center = center_X, scale = scale_X))

    if (center_X) center_X <- attr(assay(X, "transformed"), "scaled:center")
    if (scale_X) scale_X <- attr(assay(X, "transformed"), "scaled:scale")

    wgcna_res <- run_wgcna(assay(X, "transformed"),
        return_full_output = FALSE, ...
    )
    reduced_set <- .se_to_me(
        X,
        reduced = wgcna_res$E,
        loadings = wgcna_res$L,
        assignments = wgcna_res$assignments,
        center_X = center_X,
        scale_X = scale_X,
        dendrogram = wgcna_res$dendrogram,
        threshold = wgcna_res$threshold
    )

    return(reduced_set)
}

#' Creates a ModularExperiment from a SummarizedExperiment
#'
#' Helper function for transforming a
#' \link{ReducedExperiment}[ModularExperiment] into a
#' \link{SummarizedExperiment}[SummarizedExperiment]
#'
#' @param se A \link{SummarizedExperiment}[SummarizedExperiment] object.
#'
#' @param reduced Data to be passed to the `reduced` slot.
#'
#' @param loadings Data to be passed to the `loadings` slot.
#'
#' @param assignments Data to be passed to the `assignments` slot.
#'
#' @param center_X Data to be passed to the `center_X` slot.
#'
#' @param scale_X Data to be passed to the `scale_X` slot.
#'
#' @param dendrogram Data to be passed to the `dendrogram` slot.
#'
#' @param threshold Data to be passed to the `threshold` slot.
#'
#' @noRd
#' @keywords internal
.se_to_me <- function(
        se,
        reduced,
        loadings,
        assignments,
        center_X,
        scale_X,
        dendrogram = NULL,
        threshold = NULL) {
    return(ModularExperiment(
        reduced = reduced, loadings = loadings, assignments = assignments,
        center = center_X, scale = scale_X,
        dendrogram = dendrogram, threshold = threshold,
        assays = assays(se), rowData = rowData(se),
        colData = colData(se), metadata = S4Vectors::metadata(se)
    ))
}

#' Run WGCNA for a data matrix
#'
#' Runs WGCNA. Largely a wrapper for the \link[WGCNA]{blockwiseModules}
#' function. Additionally applies \link[WGCNA]{pickSoftThreshold} to
#' aid in the selection of the soft thresholding power, reformats data
#' into a format convenientg for creating a
#' \link[ReducedExperiment]{ModularExperiment} object and changes module names
#' from colours to numbers (default).
#'
#' The function also stores the loadings matrices generated when PCA is
#' performed for each module to calculate eigengenes. These loadings can be
#' used to quantify the alignment of genes with the module eigengene,
#' can be used to recalculate the reduced data matrix (eigengenes) and these
#' loadings can be used to transform new datasets.
#'
#' @param X A matrix with features as rows and columns as samples.
#'
#' @param powers The soft-thresholding power(s) to test, see
#' \link[WGCNA]{pickSoftThreshold} for more details.
#'
#' @param min_r_squared The correlation threshold for selection of the
#' soft-thresholding power, see \link[WGCNA]{pickSoftThreshold} for more
#' details.
#'
#' @param max_mean_connectivity A threshold for the maximal mean connectivity,
#' used when selecting the soft-thresholding power, see
#' \link[WGCNA]{pickSoftThreshold} for more details.
#'
#' @param corType The type of correlation to be used to generate a correlation
#' matrix during network formation. One of "pearson" (\link[WGCNA]{cor}) and
#' "bicor" (\link[WGCNA]{bicor}).
#'
#' @param networkType The type of network to be generated, either "unsigned",
#' "signed" or "signed hybrid". See \link[WGCNA]{blockwiseModules} for
#' more details.
#'
#' @param module_labels Specifies whether the modules should be named based on
#' "numbers" or "colours. If `module_labels` is set to "numbers", then
#' "module_0" represents unclustered genes, whereas if it is set to "colours"
#' then "grey" represents unclustered genes.
#'
#' @param maxBlockSize The chunk size (in terms of the number of features/genes)
#' to process the data. See \link[WGCNA]{blockwiseModules} for
#' more details. The default (30000) should process standard transcriptomic
#' datasets in a single chunk. Results may differ if the number of features
#' exceeds the chunk size. Lower values of this parameter may use less memory
#' to calculate networks.
#'
#' @param verbose The verbosity, passed to \link[WGCNA]{blockwiseModules}.
#'
#' @param return_full_output If FALSE (default), returns the results specified
#' below in `returns`. Else, returns additional information, including
#' "run_wgcna_output") the results specified below in `returns`; "original_E")
#' the original eigengene matrix calculated by \link[WGCNA]{blockwiseModules};
#' "blockwise_modules_output") the object returned by
#' \link[WGCNA]{blockwiseModules}; and "pick_soft_threshold_output") the
#' output of \link[WGCNA]{pickSoftThreshold}.
#'
#' @param scale_reduced If TRUE, the reduced data (eigengenes) are standardised
#' to have a mean of 0 and a standard deviation of 1.
#'
#' @param ... Additional arguments to be passed to
#' \link[WGCNA]{blockwiseModules}.
#'
#' @details
#' Note that if `module_labels` is set to "numbers", then
#' "module_0" represents unclustered genes, whereas if it is set to "colours"
#' then "grey" represents unclustered genes.
#'
#' Note that if `powers` is a range, the power to use will be selected
#' automatically based on the output of \link[WGCNA]{pickSoftThreshold}. We
#' recommend that users avoid relying on the automatic threshold detection
#' and instead consider the best soft thresholding power carefully, for instance
#' using the plots described in the `WGCNA` tutorials.
#'
#' @returns Returns a list containing:
#' \describe{
#'  \item{"E"}{The reduced data (eigengenes).}
#'  \item{"L"}{The module loadings. This represents the values of the PCA
#'  rotation matrix for the first principal component generated for each
#'  module.}
#'  \item{"assignments"}{A named vector representing the assignments of
#'  genes to modules.}
#' }
#'
#' @export
run_wgcna <- function(
        X, powers = 1:30, min_r_squared = 0.85,
        max_mean_connectivity = 100, corType = "pearson",
        networkType = "signed", module_labels = "numbers", maxBlockSize = 30000,
        verbose = 0, return_full_output = FALSE, scale_reduced = TRUE, ...) {
    .max_block_size_check(maxBlockSize, nrow(X))
    cor <- corFnc <- .get_cor_fn(corType) # Get correlation function

    threshold <- .select_threshold(t(X),
        min_r_squared = min_r_squared, maxBlockSize = maxBlockSize,
        verbose = verbose, max_mean_connectivity = max_mean_connectivity,
        powers = powers, corFnc = corFnc, networkType = networkType
    ) # Suggest best threshold

    power <- unique(threshold$fitIndices$selected_power)

    bwms <- WGCNA::blockwiseModules(
        t(X),
        power = power, corType = corType, networkType = networkType,
        maxBlockSize = maxBlockSize, verbose = verbose, ...
    ) # Apply WGCNA pipeline

    d <- if (length(bwms$dendrograms) == 1) bwms$dendrograms[[1]] else NULL
    wgcna_res <- list(
        "assignments" = bwms$colors, "E" = bwms$MEs,
        "dendrogram" = d, "threshold" = threshold$fitIndices
    )

    original_E <- wgcna_res$E
    colnames(wgcna_res$E) <- gsub("ME", "", colnames(wgcna_res$E))

    if (module_labels == "numbers") { # Convert colours to numbers
        converter <- .colors2numbers(wgcna_res$assignments)
        wgcna_res$assignments <- vapply(
            wgcna_res$assignments, converter,
            FUN.VALUE = 1
        )
        colnames(wgcna_res$E) <- vapply(
            colnames(wgcna_res$E), converter,
            FUN.VALUE = 1
        )
    } else if (!module_labels %in% c("colors", "colours")) {
        stop("Value of `module_labels` does not correspond to a valid option")
    }

    wgcna_res$assignments <- stats::setNames(
        names(wgcna_res$assignments),
        paste0("module_", wgcna_res$assignments)
    )
    colnames(wgcna_res$E) <- paste0("module_", colnames(wgcna_res$E))
    wgcna_res$E <- wgcna_res$E[, order(colnames(wgcna_res$E))]

    original_E <- wgcna_res$E
    recalculated_E <- .calculate_eigengenes(
        X,
        colnames(wgcna_res$E),
        wgcna_res$assignments,
        realign = TRUE
    )

    wgcna_res$E <- recalculated_E$reduced
    wgcna_res$L <- recalculated_E$loadings
    if (scale_reduced) wgcna_res$E <- scale(wgcna_res$E)

    if (return_full_output) {
        return(list(
            "run_wgcna_output" = wgcna_res,
            "original_E" = original_E,
            "blockwise_modules_output" = bwms,
            "pick_soft_threshold_output" = threshold
        ))
    } else {
        return(wgcna_res)
    }
}

#' Based on a string, determine the WGCNA correlation function to use
#'
#' @noRd
#' @keywords internal
.get_cor_fn <- function(corType) {
    if (corType == "pearson") {
        return(WGCNA::cor)
    } else if (corType == "bicor") {
        return(WGCNA::bicor)
    } else {
        stop("`corType` must be one of 'pearson', 'bicor'")
    }
}

#' Raise warning if blockwise module detection is taking place
#'
#' @noRd
#' @keywords internal
.max_block_size_check <- function(max_block_size, n_rows) {
    if (max_block_size < n_rows) {
        warning(
            "maxBlockSize < total features, module detection will be ",
            "performed in a block-wise manner"
        )
    }
}

#' Wrapper around WGCNA::pickSoftThreshold
#'
#' @noRd
#' @keywords internal
.select_threshold <- function(X,
    min_r_squared,
    max_mean_connectivity,
    powers,
    corFnc,
    networkType,
    maxBlockSize,
    verbose) {
    threshold <- WGCNA::pickSoftThreshold(
        t(X),
        RsquaredCut = min_r_squared, powerVector = powers, corFnc = corFnc,
        networkType = networkType, blockSize = maxBlockSize, verbose = verbose
    )

    if (length(powers) > 1) {
        if (is.null(max_mean_connectivity)) {
            power <- threshold$fitIndices$powerEstimate
        } else {
            which_power <- which(
                threshold$fitIndices$SFT.R.sq > min_r_squared &
                    threshold$fitIndices$mean.k. < max_mean_connectivity
            )

            if (length(which_power) == 0) {
                stop(
                    "No power with r_squared > ", min_r_squared,
                    " and mean connectivity < ", max_mean_connectivity
                )
            }

            power <- min(threshold$fitIndices$Power[which_power])
        }
    } else if (length(powers) == 1) {
        power <- powers
    }

    threshold$fitIndices$selected_power <- power

    return(threshold)
}

#' Convert WGCNA colours to numbers
#'
#' @noRd
#' @keywords internal
.colors2numbers <- function(colors) {
    color_table <- table(colors)
    color_table <- color_table[order(color_table, decreasing = TRUE)]
    color_table <- color_table[which(names(color_table) != "grey")]

    color_table <- stats::setNames(
        seq_len(length(color_table)),
        names(color_table)
    )
    color_table <- c(color_table, "grey" = 0)

    return(function(x) {
        color_table[which(names(color_table) == x)]
    })
}

#' Calculates eigengenes for new data based on module assignments
#'
#' @noRd
#' @keywords internal
.calculate_eigengenes <- function(
        newdata,
        module_names,
        module_assignments,
        realign = TRUE) {
    red <- data.frame(row.names = colnames(newdata))
    lod <- c()

    for (m in module_names) {
        module_features <- module_assignments[names(module_assignments) == m]
        module_data <- newdata[rownames(newdata) %in% module_features, ]

        prcomp_res <- stats::prcomp(t(module_data),
            center = FALSE, scale. = FALSE, rank = 1
        )

        # Principal components may be anticorrelated, in which case change sign
        align_sign <- ifelse(
            realign,
            sign(mean(stats::cor(prcomp_res$x, t(module_data)))),
            1
        )

        red[[m]] <- prcomp_res$x * align_sign

        stopifnot(all(t(module_data) %*% prcomp_res$rotation * align_sign
        == prcomp_res$x * align_sign))

        lod <- c(lod, stats::setNames(
            prcomp_res$rotation * align_sign,
            rownames(module_data)
        ))
    }

    lod <- lod[match(rownames(newdata), names(lod))]

    return(list("reduced" = as.matrix(red), "loadings" = lod))
}

#' Calculates eigengenes for new data based on the original PCA rotation matrix
#'
#' @noRd
#' @keywords internal
.project_eigengenes <- function(
        newdata,
        module_names,
        module_assignments,
        lod,
        min_module_genes) {
    red <- data.frame(row.names = colnames(newdata))

    for (m in module_names) {
        module_genes <- module_assignments[names(module_assignments) == m]

        if (length(module_genes) < min_module_genes) next
        if (!any(rownames(newdata) %in% module_genes)) next


        module_data <- newdata[rownames(newdata) %in% module_genes, ]
        module_lod <- lod[names(lod) %in% module_genes]
        stopifnot(all(rownames(module_data) == names(module_lod)))

        red[[m]] <- t(module_data) %*% module_lod
    }

    return(as.matrix(red))
}
