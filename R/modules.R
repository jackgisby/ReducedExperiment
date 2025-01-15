#' Apply dimensionality reduction using Weighted Gene Correlation Network Analysis
#'
#' Performs Weighted gene correlation network analysis (WGCNA) and packages
#' both the input data and subsequent results into a
#' \link[ReducedExperiment]{ModularExperiment}. Calls
#' \link[ReducedExperiment]{runWGCNA} to perform the analysis; see its
#' documentation page for more information on the ICA method, parameters
#' and outputs.
#'
#' @param X Either a \link[SummarizedExperiment]{SummarizedExperiment} object
#' or a matrix containing data to be subject to WGCNA. `X` should have rows as
#' features and columns as samples.
#'
#' @param power An integer representing the soft-thresholding power to be
#' used to define modules. See the
#' \link[ReducedExperiment]{assessSoftThreshold} function for aid in
#' determining this parameter.
#'
#' @param center_X If `TRUE`, `X` is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to WGCNA.
#'
#' @param scale_X If `TRUE`, `X` is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before WGCNA.
#'
#' @param assay_name If `X` is a
#' \link[SummarizedExperiment]{SummarizedExperiment}, then this should be the
#' name of the assay to be subject to WGCNA.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{runWGCNA}.
#'
#' @returns A \link[ReducedExperiment]{ModularExperiment} is returned
#' containing the input data (i.e., the original data matrix in addition to
#' other slots if a \link[SummarizedExperiment]{SummarizedExperiment} was used
#' as input). Additionally contains the results of module analysis, stored in
#' the `reduced` and `assignments` slots. The `center_X`, `scale_X`,
#' `loadings`, `threshold` and `dendrogram` slots may also be filled depending
#' on the arguments given to `identifyModules`.
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get the airway data as a SummarizedExperiment (with a subset of features)
#' set.seed(2)
#' airway_se <- ReducedExperiment:::.getAirwayData(n_features = 500)
#'
#' # Select soft-thresholding power to use (use capture.output to hide WGCNA's prints)
#' WGCNA::disableWGCNAThreads()
#' invisible(capture.output(fit_indices <- assessSoftThreshold(airway_se)))
#' estimated_power <- fit_indices$Power[fit_indices$estimated_power]
#'
#' # Identify modules using WGCNA
#' airway_me <- identifyModules(airway_se, verbose = 0, power = estimated_power)
#' airway_me
#'
#' @seealso [ReducedExperiment::runWGCNA()],
#'     [WGCNA::blockwiseModules()],
#'     [WGCNA::pickSoftThreshold()]
#'
#' @export
identifyModules <- function(
    X, power, center_X = TRUE, scale_X = TRUE, assay_name = "normal",
    ...
) {
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

    wgcna_res <- runWGCNA(assay(X, "transformed"), power = power, ...)
    reduced_set <- .seToMe(
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
#' \link[ReducedExperiment]{ModularExperiment} into a
#' \link[SummarizedExperiment]{SummarizedExperiment}
#'
#' @param se A \link[SummarizedExperiment]{SummarizedExperiment} object.
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
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.seToMe <- function(
    se, reduced, loadings, assignments, center_X, scale_X,
    dendrogram = NULL, threshold = NULL
) {
    return(ModularExperiment(
        reduced = reduced, loadings = loadings, assignments = assignments,
        center = center_X, scale = scale_X,
        dendrogram = dendrogram, threshold = threshold,
        assays = assays(se), rowData = rowData(se),
        colData = colData(se), metadata = S4Vectors::metadata(se)
    ))
}

#' Assess soft thresholding power for WGCNA
#'
#' A wrapper around \link[WGCNA]{pickSoftThreshold}, allowing assessment and
#' automatic selection of soft-thresholding power. Extends the function to
#' accept a \link[SummarizedExperiment]{SummarizedExperiment} as input and
#' additionally considers mean connectivity when selecting the
#' soft-thresholding power to recommend.
#'
#' @param max_mean_connectivity The maximal mean connectivity required. Used
#' to select the soft-thresholding power.
#'
#' @param cor_type The type of correlation to be used to generate a correlation
#' matrix during network formation. One of "pearson" (\link[WGCNA]{cor}) and
#' "bicor" (\link[WGCNA]{bicor}).
#'
#' @param maxBlockSize The chunk size (in terms of the number of features/genes)
#' to process the data. The default (30000) should process standard
#' transcriptomic
#' datasets in a single chunk. Results may differ if the number of features
#' exceeds the chunk size. Lower values of this parameter will use less memory
#' to calculate networks.
#'
#' @param ... Additional arguments to be passed to
#' \link[WGCNA]{pickSoftThreshold}.
#'
#' @inheritParams WGCNA::pickSoftThreshold
#' @inheritParams identifyModules
#'
#' @returns Returns a `data.frame`,
#' generated by \link[WGCNA]{pickSoftThreshold},
#' with scale free topology fitting indices and connectivity statistics.
#' Additionally contains a column, `estimated_power`, indicating the
#' recommended power to use (see `details`). We suggest manually considering
#' suitability of the soft-thresholding power rather than solely relying on this
#' automated approach.
#'
#' @details
#' The \link[WGCNA]{pickSoftThreshold} function estimates the power by
#' selecting the lowest value with
#' a minimum scale free topology fitting index exceeding `RsquaredCut`.
#' The `assessSoftThreshold` function mirrors this behaviour when
#' `max_mean_connectivity` is `NULL`. When `max_mean_connectivity` is
#' specified, however, we additionally require that the selected power
#' does not exceed this connectivity threshold.
#'
#' @examples
#' # Get the airway data as a SummarizedExperiment (with a subset of features)
#' set.seed(2)
#' airway_se <- ReducedExperiment:::.getAirwayData(n_features = 500)
#'
#' # Select soft-thresholding power to use (use capture.output to hide WGCNA's prints)
#' WGCNA::disableWGCNAThreads()
#' invisible(capture.output(fit_indices <- assessSoftThreshold(airway_se)))
#'
#' print(fit_indices)
#' print(paste0("Estimated power: ", fit_indices$Power[fit_indices$estimated_power]))
#'
#' @seealso [WGCNA::pickSoftThreshold()],
#'     [ReducedExperiment::runWGCNA()]
#'
#' @author Jack Gisby
#'
#' @importFrom utils capture.output
#' @export
assessSoftThreshold <- function(
    X, assay_name = "normal", powerVector = 1:30, RsquaredCut = 0.85,
    max_mean_connectivity = 100, cor_type = "pearson", networkType = "signed",
    maxBlockSize = 30000, verbose = 0, ...
) {
    .maxBlockSizeCheck(maxBlockSize, nrow(X))
    cor <- corFnc <- .getCorFn(cor_type)  # Get correlation function

    if (inherits(X, "SummarizedExperiment")){
        X <- assay(X, assay_name)
    }

    args <- list(
        data = t(X),
        powerVector = powerVector,
        RsquaredCut = RsquaredCut,
        corFnc = corFnc,
        networkType = networkType,
        blockSize = maxBlockSize,
        verbose = verbose,
        ...
    )

    # Apply soft thresholding function (if verbose == 0 suppress the print)
    if (verbose == 0) {
        capture.output(
            threshold_output <- do.call(WGCNA::pickSoftThreshold, args),
            file = NULL
        )
    } else {
        threshold_output <- do.call(WGCNA::pickSoftThreshold, args)
    }

    # Get the output
    wgcna_power_estimate <- threshold_output$powerEstimate
    fit_indices <- threshold_output$fitIndices

    # Automatically identify the best power
    if (is.null(max_mean_connectivity)) {
        best_power <- wgcna_power_estimate
    } else {
        which_power <- which(
            fit_indices$SFT.R.sq > RsquaredCut &
            fit_indices$mean.k. < max_mean_connectivity
        )

        if (length(which_power) == 0) {
            warning(
                "No power with r_squared > ", RsquaredCut,
                " and mean connectivity < ", max_mean_connectivity
            )
            best_power <- NULL
        } else {
            best_power <- min(fit_indices$Power[which_power])
        }
    }

    if (!is.null(best_power))
        fit_indices$estimated_power <- fit_indices$Power == best_power

    return(fit_indices)
}

#' Run WGCNA for a data matrix
#'
#' Runs WGCNA. Largely a wrapper for the \link[WGCNA]{blockwiseModules}
#' function that reformats data
#' into a format convenient for creating a
#' \link[ReducedExperiment]{ModularExperiment} object and changes module names
#' from colours to numbers by default.
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
#' @param standardise_reduced If `TRUE`, the reduced data (eigengenes) are
#' standardised to have a mean of 0 and a standard deviation of 1.
#'
#' @param ... Additional arguments to be passed to
#' \link[WGCNA]{blockwiseModules}.
#'
#' @inheritParams WGCNA::blockwiseModules
#' @inheritParams identifyModules
#' @inheritParams assessSoftThreshold
#'
#' @details
#' Note that if `module_labels` is set to "numbers", then
#' "module_0" represents unclustered genes, whereas if it is set to "colours"
#' then "grey" represents unclustered genes.
#'
#' The function also stores the loadings matrices generated when PCA is
#' performed for each module to calculate eigengenes. These loadings can be
#' used to recalculate the reduced data matrix (eigengenes).
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
#' @seealso [WGCNA::blockwiseModules()],
#'     [ReducedExperiment::assessSoftThreshold()],
#'     [WGCNA::pickSoftThreshold()],
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get the airway data as a SummarizedExperiment (with a subset of features)
#' set.seed(2)
#' airway_se <- ReducedExperiment:::.getAirwayData(n_features = 500)
#'
#' # Choose an appropriate soft-thresholding power
#' WGCNA::disableWGCNAThreads()
#' fit_indices <- assessSoftThreshold(airway_se)
#' estimated_power <- fit_indices$Power[fit_indices$estimated_power]
#'
#' # Identify modules using the airway expression matrix
#' wgcna_res <- runWGCNA(
#'     assay(airway_se, "normal"),
#'     verbose = 0,
#'     power = estimated_power
#' )
#'
#' # We find just one module for this small dataset (module_0 indicates unclustered genes)
#' table(names(wgcna_res$assignments))
#'
#' @export
runWGCNA <- function(
    X, power,  cor_type = "pearson",
    networkType = "signed", module_labels = "numbers", maxBlockSize = 30000,
    verbose = 0, standardise_reduced = TRUE,
    ...
) {
    .maxBlockSizeCheck(maxBlockSize, nrow(X))
    cor <- corFnc <- .getCorFn(cor_type) # Get correlation function

    bwms <- WGCNA::blockwiseModules(
        t(X),
        power = power, cor_type = cor_type, networkType = networkType,
        maxBlockSize = maxBlockSize, verbose = verbose, ...
    ) # Apply WGCNA pipeline

    d <- if (length(bwms$dendrograms) == 1) bwms$dendrograms[[1]] else NULL
    wgcna_res <- list(
        "assignments" = bwms$colors, "E" = bwms$MEs, "dendrogram" = d
    )

    original_E <- wgcna_res$E
    colnames(wgcna_res$E) <- gsub("ME", "", colnames(wgcna_res$E))

    if (module_labels == "numbers") {  # Convert colours to numbers
        converter <- .colorsToNumbers(wgcna_res$assignments)
        wgcna_res$assignments <- vapply(
            wgcna_res$assignments, converter, FUN.VALUE = 1
        )
        colnames(wgcna_res$E) <- vapply(
            colnames(wgcna_res$E), converter, FUN.VALUE = 1
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
    recalculated_E <- .calculateEigengenes(
        X, colnames(wgcna_res$E), wgcna_res$assignments, realign = TRUE
    )

    wgcna_res$E <- recalculated_E$reduced
    wgcna_res$L <- recalculated_E$loadings
    if (standardise_reduced) wgcna_res$E <- scale(wgcna_res$E)
    return(wgcna_res)
}

#' Based on a string, determine the WGCNA correlation function to use
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.getCorFn <- function(cor_type) {
    if (cor_type == "pearson") {
        return(WGCNA::cor)
    } else if (cor_type == "bicor") {
        return(WGCNA::bicor)
    } else {
        stop("`cor_type` must be one of 'pearson', 'bicor'")
    }
}

#' Raise warning if blockwise module detection is taking place
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.maxBlockSizeCheck <- function(max_block_size, n_rows) {
    if (max_block_size < n_rows) {
        warning(
            "maxBlockSize < total features, module detection will be ",
            "performed in a block-wise manner"
        )
    }
}

#' Convert WGCNA colours to numbers
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.colorsToNumbers <- function(colors) {
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
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.calculateEigengenes <- function(
    newdata,
    module_names,
    module_assignments,
    realign = TRUE
) {
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
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.projectEigengenes <- function(
    newdata, module_names, module_assignments, lod, min_module_genes
) {
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
