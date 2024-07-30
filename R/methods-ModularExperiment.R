#' ModularExperiment: A container for the results of module analysis
#'
#' @description
#' A container inheriting from \link[ReducedExperiment]{ReducedExperiment}, that
#' contains one or more data matrices, to which module analysis has been applied
#' to identify a reduced set of features.
#'
#' @param reduced A data matrix, produced by module analysis, with rows
#' representing samples and columns representing module expression profiles.
#' Typically, this matrix contains "eigengenes" produced by the Weighted Gene
#' Correlation Network Analysis Approach, as is applied by
#' \link[ReducedExperiment]{identify_modules}.
#'
#' @param scale Either a boolean, representing whether or not the original data
#' has been scaled to unit variance, or a numeric vector indicating the
#' standard deviations of the original features (as produced by
#' \link[base]{scale}.)
#'
#' @param center Either a boolean, representing whether or not the original data
#' has been centered to have a mean of 0, or a numeric vector indicating the
#' means of the original features (as produced by
#' \link[base]{scale}.)
#'
#' @param assignments A vector of features, named according to the module to
#' which the feature belongs.
#'
#' @param loadings A numeric vector representing the loadings used to generate
#' module expression profiles. Typically, these values are obtained from the
#' rotation matrix produced by \link[stats]{prcomp}, which is used to identify
#' the first principal component of each module. The vector names represent
#' features.
#'
#' @param dendrogram Either NULL, or the dendrogram used to identify modules
#' from the original data.
#'
#' @param threshold Either NULL, or a matrix produced by
#' \link[WGCNA]{pickSoftThreshold} indicating the parameters used for network
#' construction.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{ReducedExperiment}.
#'
#' @import SummarizedExperiment
#'
#' @rdname modular_experiment
#' @export
ModularExperiment <- function(
        reduced = new("matrix"),
        scale = TRUE,
        center = TRUE,
        loadings = numeric(),
        assignments = character(),
        dendrogram = NULL,
        threshold = NULL,
        ...) {
    re <- ReducedExperiment(
        reduced = reduced,
        scale = scale,
        center = center,
        ...
    )

    return(.ModularExperiment(
        re,
        loadings = loadings,
        assignments = assignments,
        dendrogram = dendrogram,
        threshold = threshold
    ))
}

S4Vectors::setValidity2("ModularExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Assignments
    if (obj_dims[1] != length(assignments(object))) {
        msg <- c(msg, "Assignments have invalid length")
    }

    if (!all(assignments(object) == rownames(object))) {
        msg <- c(msg, "Assignments have incompatible names (rownames)")
    }

    # Loadings
    if (obj_dims[1] != length(loadings(object))) {
        msg <- c(msg, "Loadings have invalid length")
    }

    if (!all(names(loadings(object)) == rownames(object))) {
        msg <- c(msg, "Loadings have incompatible names (rownames)")
    }

    return(if (is.null(msg)) TRUE else msg)
})

#' Get and set module feature assignments
#'
#' Retrieves a vector of features (usually genes) named by the modules
#' they belong to.
#'
#' @param object \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param as_list If `TRUE`, the results are returned as a list, with an entry
#' for each module containing a list of features.
#'
#' @param value New value to replace existing assignments.
#'
#' @rdname module_assignments
#' @name assignments
#' @aliases assignments<-
#' @export assignments
NULL

#' @rdname module_assignments
#' @export
setMethod("assignments", "ModularExperiment", function(object, as_list = FALSE) {
    if (as_list) {
        a <- list()
        for (comp in componentNames(object)) {
            a[[comp]] <-
                assignments(object)[which(names(assignments(object)) == comp)]
        }
    } else {
        a <- object@assignments
    }

    return(a)
})

#' @rdname module_assignments
#' @export
setReplaceMethod("assignments", "ModularExperiment", function(object, value) {
    object@assignments <- value
    validObject(object)
    return(object)
})

#' Get and set loadings
#'
#' Method for getting and setting loadings for a
#' \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object or an
#' object that inherits from this class.
#'
#' @param scale_loadings If `TRUE`, loadings will be scaled to have a standard
#' deviation of 0. If the loadings are a matrix, this operation is performed
#' column-wise.
#'
#' @param center_loadings If `TRUE`, loadings will be centered to have a mean
#' of 0. If the loadings are a matrix, this operation is performed
#' column-wise.
#'
#' @param abs_loadings If `TRUE`, the absolute values of the loadings will be
#' returned.
#'
#' @param value New value to replace existing loadings.
#'
#' @details
#' #' If `object` is a \link[ReducedExperiment]{FactorisedExperiment}, the
#' loadings matrix will be returned, with features as rows and reduced
#' components as columns.
#'
#' If `object` is a \link[ReducedExperiment]{ModularExperiment}, the loadings
#' will be returned as a vector, with a value for each feature (usually genes).
#'
#' When available, the module loadings provide the values of the rotation matrix
#' (usually generated by \link[stats]{prcomp}) used to calculate the
#' sample-level module vectors available in the `reduced` slot. Normally, these
#' loadings are calculated for each module separately, so their values are
#' not comparable across modules.
#'
#' @rdname loadings
#' @name loadings
#' @aliases loadings<-
#' @export loadings
NULL

#' @rdname loadings
#' @export
setMethod("loadings", "ModularExperiment", function(
        object,
        scale_loadings = FALSE,
        center_loadings = FALSE,
        abs_loadings = FALSE) {
    l <- scale(
        object@loadings,
        scale = scale_loadings,
        center = center_loadings
    )
    if (abs_loadings) l <- abs(l)
    return(l[, 1])
})

#' @rdname loadings
#' @export
setReplaceMethod("loadings", "ModularExperiment", function(object, value) {
    object@loadings <- value
    validObject(object)
    return(object)
})

#' @rdname feature_names
#' @export
setReplaceMethod("names", "ModularExperiment", function(x, value) {
    x@assignments <- stats::setNames(value, names(x@assignments))
    x@loadings <- stats::setNames(x@loadings, value)

    x <- callNextMethod(x, value)
    validObject(x)

    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("featureNames", "ModularExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("rownames", "ModularExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' @rdname component_names
#' @export
setReplaceMethod("componentNames", "ModularExperiment", function(object, value) {
    curr_names <- colnames(object@reduced)
    object <- callNextMethod(object, value)
    new_names <- colnames(object@reduced)

    for (i in seq_len(length(curr_names))) {
        names(object@assignments)[
            which(names(object@assignments) == curr_names[i])
        ] <- new_names[i]
    }

    validObject(object)
    return(object)
})

#' @rdname component_names
#' @export
setMethod("moduleNames", "ModularExperiment", function(object) {
    return(componentNames(object))
})

#' @rdname component_names
#' @export
setReplaceMethod("moduleNames", "ModularExperiment", function(object, value) {
    componentNames(object) <- value
    return(object)
})

#' @rdname individual_dims
#' @export
setMethod("nModules", "ModularExperiment", function(object) {
    dim(object)[3]
})

#' Get the dendrogram stored in a ModularExperiment
#'
#' @param object \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param value New value to replace existing dendrogram.
#'
#' @rdname module_dendrogram
#' @name dendrogram
#' @aliases dendrogram<-
#' @export dendrogram
NULL

#' @rdname module_dendrogram
#' @export
setMethod("dendrogram", "ModularExperiment", function(object) {
    return(object@dendrogram)
})

#' @rdname module_dendrogram
#' @export
setReplaceMethod("dendrogram", "ModularExperiment", function(object, value) {
    object@dendrogram <- value
    return(object)
})

# Required for dollarsign autocomplete of colData columns
.DollarNames.ModularExperiment <- function(x, pattern = "") {
    grep(pattern, colnames(colData(x)), value = TRUE)
}

#' @rdname slice
#' @export
setMethod(
    "[", c("ModularExperiment", "ANY", "ANY", "ANY"),
    function(x, i, j, k, ..., drop = FALSE) {
        object <- x

        if (1L != length(drop) || (!missing(drop) && drop)) {
            warning("'drop' ignored '[,", class(object), ",ANY,ANY-method'")
        }

        assignments <- object@assignments
        lod <- object@loadings

        if (!missing(i)) {
            if (is.character(i)) {
                fmt <- paste0(
                    "<", class(object),
                    ">[i,] index out of bounds: %s"
                )
                i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    i, rownames(object), fmt
                )
            }

            i <- as.vector(i)
            assignments <- assignments[i, drop = FALSE]
            lod <- lod[i, drop = FALSE]
        }

        out <- callNextMethod(object, i, j, k, ...)
        BiocGenerics:::replaceSlots(
            out,
            loadings = lod,
            assignments = assignments,
            check = FALSE
        )
    }
)

# Same features, different samples
#' @rdname cbind
#' @export
setMethod("cbind", "ModularExperiment", function(..., deparse.level = 1) {
    args <- list(...)

    loadings_assignments_equal <- vapply(args, function(re) {
        return(identical(re@loadings, args[[1]]@loadings)
               & identical(re@assignments, args[[1]]@assignments))
    },
        FUN.VALUE = FALSE
    )

    if (!all(loadings_assignments_equal)) {
        stop("Row bind expects loadings and assignments slots are equal")
    }

    args[["deparse.level"]] <- deparse.level

    return(do.call(callNextMethod, args))
})

#' @rdname enrichment
#' @export
setMethod("runEnrich", c("ModularExperiment"), function(
        object,
        method = "overrepresentation",
        feature_id_col = "rownames",
        as_dataframe = FALSE,
        ...) {
    if (method == "overrepresentation") {
        if (feature_id_col != "rownames") {
            names(object) <-
                rowData(object)[[feature_id_col]]
        }

        modules <- assignments(object, as_list = TRUE)

        enrich_res <- reduced_oa(modules, ...)
    } else {
        stop("Enrichment method not recognised")
    }

    if (as_dataframe) {
        enrich_res <- lapply(enrich_res, function(object) {
            object@result
        })
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})

#' Plot a dendrogram stored in a ModularExperiment
#'
#' Gets the dendrogram and plots it using
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @param object \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param groupLabels Module label axis label. See
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @param dendroLabels If TRUE, shows feature names in the dendrogram. See
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @param hang The fraction of the plot height by which labels should hang
#' below the rest of the plot. See \link[stats]{plot.hclust}.
#'
#' @param addGuide If TRUE, adds vertical guide lines to the dendrogram. See
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @param guideHang The fraction of the dendrogram's height to leave between
#' the top end of the guide line and the dendrogram merge height. See
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @param color_func Function for converting module names to colors. Only used
#' if `modules_are_colors` is FALSE
#'
#' @param modules_are_colors If TRUE, expects the module names to be colors.
#' Else, assumes that module names are are numbers that can be converted into
#' colours by `color_func`.
#'
#' @param ... Additional arguments to be passed to
#' \link[WGCNA]{plotDendroAndColors}.
#'
#' @rdname plotDendro
#' @name plotDendro
#' @export plotDendro
NULL

#' @rdname plotDendro
#' @export
setMethod(
    "plotDendro", c("ModularExperiment"),
    function(
        object, groupLabels = "Module colors", dendroLabels = FALSE,
        hang = 0.03, addGuide = TRUE, guideHang = 0.05,
        color_func = WGCNA::labels2colors, modules_are_colors = FALSE, ...) {
        if (!modules_are_colors) {
            colors <- as.numeric(gsub(
                "module_", "",
                names(assignments(object))
            ))
            colors <- color_func(colors)
        }

        WGCNA::plotDendroAndColors(dendrogram(object), colors,
            groupLabels = groupLabels,
            dendroLabels = dendroLabels, hang = hang,
            addGuide = addGuide, guideHang = guideHang,
            ...
        )
    }
)

#' Project new data using pre-defined factors
#'
#' @description
#' Calculates eigengenes for modules in new data. If in `project` mode,
#' functions in a similar fashion to the `predict` method of
#' \link[stats]{prcomp}. Else, eigengenes are calculated from scratch using PCA,
#' in a similar manner to the \link[WGCNA]{moduleEigengenes} function.
#'
#' @param object A \link[ReducedExperiment]{ModularExperiment} object. The
#' `loadings` slot of this class will be used for projection. Additionally,
#' by default, the `scale` and `center` slots are used to apply the original
#' transformation to the new data.
#'
#' @param newdata New data for eigengenes to be calculates in. Must be a
#' `data.frame` or `matrix` with features as rows and samples as columns, or a
#' \link[SummarizedExperiment]{SummarizedExperiment} object. Assumes that the
#' rows of `newdata` match those of the
#' \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param project Whether to perform projection (i.e., using PCA rotation matrix
#' from the original data to calculate modules) or calculate eigengenes from
#' scratch in the new data (i.e., performing PCA for each module in `newdata`).
#'
#' @param scale_reduced Whether or not the reduced data should be scaled
#' after calculation.
#'
#' @param scale_newdata Controls whether the `newdata` are scaled. If NULL,
#' performs scaling based on the \link[ReducedExperiment]{ModularExperiment}
#' object's `scale` slot. The value of this argument will be passed to the
#' `scale` argument of \link[base]{scale}.
#'
#' @param center_newdata Controls whether the `newdata` are centered If NULL,
#' performs centering based on the \link[ReducedExperiment]{ModularExperiment}
#' object's `center` slot. The value of this argument will be passed to the
#' `center` argument of \link[base]{scale}.
#'
#' @param assay_name If a \link[SummarizedExperiment]{SummarizedExperiment}
#' object is passed as new data, this argument indicates which assay should be
#' used for projection.
#'
#' @param realign If `project` is TRUE, this argument is ignored. Else, controls
#' whether eigengenes are realigned after PCA is performed to ensure the
#' resultant signatures are positively correlated with average expression of the
#' module. Similar to the `align` argument of \link[WGCNA]{moduleEigengenes}.
#'
#' @param min_module_genes If `project` is FALSE, this argument is ignores.
#' Else, controls the minimum number of genes required in a module for
#' projection. Projected eigengenes are not calculated for modules with sizes
#' below this threshold.
#'
#' @param return_loadings If True, additionally returns the feature loadings for
#' the eigengenes.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{calcEigengenes}.
#'
#' @returns If return_loadings is True, returns a list with the "reduced" matrix
#' and "loadings" vector (one value per feature). If False, returns only the
#' reduced matrix.
#'
#' The reduced matrix has samples as rows and modules as columns. If
#' `newdata` was a `matrix` or `data.frame`, this will be returned as a matrix.
#' If a \link[SummarizedExperiment]{SummarizedExperiment} object was passed
#' instead, then a If a \link[ReducedExperiment]{ModularExperiment}
#' object will be created containing this matrix in its `reduced` slot.
#'
#' @seealso \code{\link[ReducedExperiment]{projectData}},
#' \code{\link[WGCNA]{moduleEigengenes}}
#'
#' @rdname calcEigengenes
#' @name calcEigengenes
#' @export calcEigengenes
NULL

#' @rdname calcEigengenes
#' @export
setMethod("calcEigengenes", c("ModularExperiment", "matrix"), function(
        object,
        newdata,
        project = TRUE,
        scale_reduced = TRUE,
        return_loadings = FALSE,
        scale_newdata = NULL,
        center_newdata = NULL,
        realign = TRUE,
        min_module_genes = 10) {
    if (!identical(rownames(object), rownames(newdata))) {
        stop("Rownames of x do not match those of newdata")
    }

    # apply known vectors for scaling and centering (returned as attributes
    # by `scale`)
    if (is.null(scale_newdata)) scale_newdata <- object@scale
    if (is.null(center_newdata)) center_newdata <- object@center

    newdata <- t(scale(
        t(newdata),
        scale = scale_newdata,
        center = center_newdata
    ))

    if (project) {
        red <- .project_eigengenes(
            newdata,
            moduleNames(object),
            assignments(object),
            loadings(object),
            min_module_genes = min_module_genes
        )
        eigengenes <- list("reduced" = as.matrix(red), "
                           loadings" = loadings(object))
    } else {
        eigengenes <- .calculate_eigengenes(
            newdata,
            moduleNames(object),
            assignments(object),
            realign = realign
        )
    }

    if (scale_reduced) eigengenes$red <- scale(eigengenes$red)

    if (return_loadings) {
        return(eigengenes)
    } else {
        return(eigengenes$red)
    }
})

#' @rdname calcEigengenes
#' @export
setMethod("calcEigengenes", c("ModularExperiment", "data.frame"), function(
        object,
        newdata,
        project = TRUE,
        scale_reduced = TRUE,
        return_loadings = FALSE,
        scale_newdata = NULL,
        center_newdata = NULL,
        realign = TRUE,
        min_module_genes = 10) {
    return(calcEigengenes(object, as.matrix(newdata),
        project = project, return_loadings = return_loadings,
        scale_newdata = scale_newdata, center_newdata = center_newdata,
        realign = realign, scale_reduced = scale_reduced,
        min_module_genes = min_module_genes
    ))
})

#' @rdname calcEigengenes
#' @export
setMethod(
    "calcEigengenes", c("ModularExperiment", "SummarizedExperiment"),
    function(
        object,
        newdata,
        project = TRUE,
        scale_reduced = TRUE,
        assay_name = "normal",
        scale_newdata = NULL,
        center_newdata = NULL,
        realign = TRUE,
        min_module_genes = 10) {
        eig <- calcEigengenes(object, assay(newdata, assay_name),
            project = project, return_loadings = FALSE,
            scale_newdata = scale_newdata, center_newdata = center_newdata,
            realign = realign,
            scale_reduced = scale_reduced, min_module_genes = min_module_genes
        )

        return(.se_to_me(newdata,
            reduced = as.matrix(eig),
            loadings = loadings(object),
            assignments = assignments(object),
            center_X = object@center, scale_X = object@scale
        ))
    }
)

#' @rdname calcEigengenes
#' @export
setMethod("predict", c("ModularExperiment"), function(object, newdata, ...) {
    return(calcEigengenes(object, newdata, ...))
})

#' Get correlation of features with module eigengenes
#'
#' @param object \link[ReducedExperiment]{ModularExperiment} object.
#'
#' @param assay_name The name of the assay to be used for calculation of
#' module centrality.
#'
#' @param feature_id_col The column in `rowData(object)` that will be used as a
#' feature ID. Setting this to "rownames" (default) instead uses
#' `rownames(object)`.
#'
#' Results indicate correlation (r) and squared correlation (rsq) with the
#' module eigengene. Feature ranks within each module are also returned based
#' on these statistics.
#'
#' @rdname getCentrality
#' @name getCentrality
#' @export getCentrality
NULL

#' @rdname getCentrality
#' @export
setMethod("getCentrality", c("ModularExperiment"), function(
        object,
        assay_name = "normal",
        feature_id_col = "rownames") {
    # Get module membership (correlation with eigengene)
    signed_kme <- WGCNA::signedKME(
        t(assay(object, assay_name)),
        reduced(object)
    )
    colnames(signed_kme) <- componentNames(object)

    stopifnot(all(rownames(signed_kme) == rownames(object)))
    if (feature_id_col != "rownames") {
        rownames(signed_kme) <-
            rowData(object)[[feature_id_col]]
    }

    # Transform into a dataframe with relevant statistics
    module_features <- data.frame()
    for (m in componentNames(object)) {
        which_features <- which(names(assignments(object)) == m)

        module_kme <- data.frame(
            module = m,
            feature = rownames(signed_kme)[which_features],
            r = signed_kme[[m]][which_features]
        )

        module_kme$rsq <- module_kme$r**2
        module_kme$rank_r <- rank(1 - module_kme$r)
        module_kme$rank_rsq <- rank(1 - module_kme$rsq)

        module_features <- rbind(module_features, module_kme)
    }

    module_features <- module_features[order(module_features$rsq,
        decreasing = TRUE
    ), ]
    module_features <- module_features[order(module_features$module), ]

    return(module_features)
})
