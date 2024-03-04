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
#' \link[identify_modules].
#'
#' @param assignments A vector of features, named according to the module to which the
#' feature belongs.
#'
#' @param loadings A numeric vector representing the loadings used to generate
#' module expression profiles. Typically, these values are obtained from the
#' rotation matrix produced by \link[stats]{prcomp}, which is used to identify
#' the first principal component of each module. The vector names represent
#' features.
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
#' @param dendrogram Either NULL, or the dendrogram used to identify modules from the
#' original data.
#'
#' @param threshold Either NULL, or a matrix produced by
#' \link[WGCNA]{pickSoftThreshold} indicating the parameters used for network
#' construction.
#'
#' A class to represent the results of module analysis
#'
#' @import SummarizedExperiment
#'
#' @export
ModularExperiment <- function(loadings=numeric(), assignments=character(), dendrogram=NULL, threshold=NULL, ...)
{
    re <- ReducedExperiment(...)
    return(.ModularExperiment(re, loadings=loadings, assignments=assignments, dendrogram=dendrogram, threshold=threshold))
}

S4Vectors::setValidity2("ModularExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Assignments
    if (obj_dims[1] != length(assignments(object)))
        msg <- c(msg, "Assignments have invalid length")

    if (!all(assignments(object) == rownames(object)))
        msg <- c(msg, "Assignments have incompatible names (rownames)")

    # Loadings
    if (obj_dims[1] != length(loadings(object)))
        msg <- c(msg, "Loadings have invalid length")

    if (!all(names(loadings(object)) == rownames(object)))
        msg <- c(msg, "Loadings have incompatible names (rownames)")

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("assignments", "ModularExperiment", function(x, as_list=FALSE) {
    if (as_list)  {
        a <- list()
        for (comp in componentNames(x)) {
            a[[comp]] <- assignments(x)[which(names(assignments(x)) == comp)]
        }
    } else {
        a <- x@assignments
    }

    return(a)
})
setReplaceMethod("assignments", "ModularExperiment", function(x, value) {
    x@assignments <- value
    validObject(x)
    return(x)
})

setMethod("loadings", "ModularExperiment", function(x, scale_loadings=FALSE, center_loadings=FALSE, abs_loadings=FALSE) {
    l <- scale(x@loadings, scale=scale_loadings, center=center_loadings)
    if (abs_loadings) l <- abs(l)
    return(l[,1])
})
setReplaceMethod("loadings", "ModularExperiment", function(x, value) {
    x@loadings <- value
    validObject(x)
    return(x)
})

setReplaceMethod("names", "ModularExperiment", function(x, value) {
    x@assignments <- setNames(value, names(x@assignments))
    x@loadings <- setNames(x@loadings, value)
    x <- callNextMethod(x, value)
    validObject(x)
    return(x)
})
setReplaceMethod("featureNames", "ModularExperiment", function(x, value) {
    names(x) <- value
    return(x)
})
setReplaceMethod("rownames", "ModularExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

setReplaceMethod("componentNames", "ModularExperiment", function(x, value) {

    curr_names <- colnames(x@reduced)
    x <- callNextMethod(x, value)
    new_names <- colnames(x@reduced)

    for (i in 1:length(curr_names)) {
        names(x@assignments)[which(names(x@assignments) == curr_names[i])] <- new_names[i]
    }

    validObject(x)
    return(x)
})

setMethod("moduleNames", "ModularExperiment", function(x) {
    return(componentNames(x))
})

setReplaceMethod("moduleNames", "ModularExperiment", function(x, value) {
    componentNames(x) <- value
    return(x)
})

setMethod("nModules", "ModularExperiment", function(x) {dim(x)[3]})

setMethod("dendrogram", "ModularExperiment", function(x) {
    return(x@dendrogram)
})

setReplaceMethod("dendrogram", "ModularExperiment", function(x, value) {
    dendrogram(x) <- value
    return(x)
})

setMethod("[", c("ModularExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    assignments <- x@assignments
    lod <- x@loadings

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
              i, rownames(x), fmt
            )
        }

        i <- as.vector(i)
        assignments <- assignments[i, drop=FALSE]
        lod <- lod[i, drop=FALSE]
    }

    out <- callNextMethod(x, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, loadings=lod, assignments=assignments, check=FALSE)
})

# Same features, different samples
setMethod("cbind", "ModularExperiment", function(..., deparse.level=1) {
    args <- list(...)

    loadings_assignments_equal <- sapply(args, function(re) {
        return(identical(re@loadings, args[[1]]@loadings) & identical(re@assignments, args[[1]]@assignments))
    })

    if (!all(loadings_assignments_equal))
        stop("Row bind expects loadings and assignments slots are equal")

    args[["deparse.level"]] <- deparse.level

    return(do.call(callNextMethod, args))
})

setMethod("runEnrich", c("ModularExperiment"),
          function(x, method="overrepresentation", feature_id_col="rownames", as_dataframe=FALSE, ...)
{
    if (method == "overrepresentation") {

        if (feature_id_col != "rownames") names(x) <- rowData(x)[[feature_id_col]]

        modules <- assignments(x, as_list=TRUE)

        enrich_res <- reduced_oa(modules, ...)

    } else {
        stop("Enrichment method not recognised")
    }

    if (as_dataframe)  {
        enrich_res <- lapply(enrich_res, function(x) {x@result})
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})

setMethod("plotDendro", c("ModularExperiment"),
          function(x, groupLabels = "Module colors", dendroLabels = FALSE,
                       hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                       color_func = WGCNA::labels2colors) {

    colors <- gsub("module_", "", names(assignments(x)))

    is_color <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)),
                     error = function(e) FALSE)
        })
    }

    if (!is.null(color_func) & !all(is_color(colors))) {
        colors <- color_func(colors)
    }

    WGCNA::plotDendroAndColors(dendrogram(x), colors,
                           groupLabels = groupLabels,
                           dendroLabels = dendroLabels, hang = hang,
                           addGuide = addGuide, guideHang = guideHang)
})

#' Project new data using pre-defined factors
#'
#' @description
#' Calculates eigengenes for modules in new data. If in `project` mode, functions
#' in a similar fashion to the `predict` method of \link[stats]{prcomp}. Else,
#' eigengenes are calculated from scratch using PCA, in a similar manner to the
#' \link[WGCNA]{moduleEigengenes} function.
#'
#' @param x A \link[ModularExperiment] object. The `loadings` slot of
#' this class will be used for projection. Additionally, by default, the `scale`
#' and `center` slots are used to apply the original transformation to the
#' new data.
#'
#' @param newdata New data for eigengenes to be calculates in. Must be a
#' `data.frame` or `matrix` with features as rows and samples as columns, or a
#' \link[SummarizedExperiment]{SummarizedExperiment} object. Assumes that the
#' rows of `newdata` match those of the \link[ReducedExperiment]{ModularExperiment}
#' object.
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
#' whether eigengenes are realigned after PCA is performed to ensure the resultant
#' signatures are positively correlated with average expression of the module.
#' Similar to the `align` argument of \link[WGCNA]{moduleEigengenes}.
#'
#' @param min_module_genes If `project` is FALSE, this argument is ignores. Else,
#' controls the minimum number of genes required in a module for projection. Projected
#' eigengenes are not calculated for modules with sizes below this threshold.
#'
#' @returns Calculates a matrix with samples as rows and modules as columns. If
#' `newdata` was a `matrix` or `data.frame`, this will be returned as a matrix.
#' If a \link[SummarizedExperiment]{SummarizedExperiment} object was passed
#' instead, then a If a \link[SummarizedExperiment]{ModularExperiment}
#' object will be created containing this matrix in its `reduced` slot.
#'
#' @seealso [ReducedExperiment::projectData()], [WGCNA::moduleEigengenes()]
#'
#' @rdname calcEigengenes
#' @export
setMethod("calcEigengenes", c("ModularExperiment", "matrix"),
          function(x, newdata, project=TRUE, scale_reduced=TRUE, return_loadings=FALSE, scale_newdata=NULL, center_newdata=NULL, realign=TRUE, min_module_genes=10) {

    if (!identical(rownames(x), rownames(newdata)))
        stop("Rownames of x do not match those of newdata")

    # apply known vectors for scaling and centering (returned as attributes by `scale`)
    if (is.null(scale_newdata)) scale_newdata <- x@scale
    if (is.null(center_newdata)) center_newdata <- x@center

    newdata <- t(scale(t(newdata), scale=scale_newdata, center=center_newdata))

    if (project) {

        red <- .project_eigengenes(newdata, moduleNames(x), assignments(x), loadings(x), min_module_genes=min_module_genes)
        eigengenes <- list("reduced" = as.matrix(red), "loadings" = loadings(x))

    } else {

        eigengenes <- .calculate_eigengenes(newdata, moduleNames(x), assignments(x), realign=realign)

        # eig <- WGCNA::moduleEigengenes(t(newdata), setNames(names(assignments(x)), assignments(x)), ...)
        # colnames(eig$eigengenes) <- gsub("ME", "", colnames(eig$eigengenes))
        # red <- eig$eigengenes
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
        x, newdata, project=TRUE, scale_reduced=TRUE, return_loadings=FALSE, scale_newdata=NULL, center_newdata=NULL, realign=TRUE, min_module_genes=10) {

    return(calcEigengenes(x, as.matrix(newdata), project=project, return_loadings=return_loadings,
                          scale_newdata=scale_newdata, center_newdata=center_newdata, realign=realign,
                          scale_reduced=scale_reduced, min_module_genes=min_module_genes))
})

#' @rdname calcEigengenes
#' @export
setMethod("calcEigengenes", c("ModularExperiment", "SummarizedExperiment"),
          function(x, newdata, project=TRUE, scale_reduced=TRUE, assay_name="normal", scale_newdata=NULL, center_newdata=NULL, realign=TRUE, min_module_genes=10) {

    eig <- calcEigengenes(x, assay(newdata, assay_name), project=project, return_loadings=FALSE,
                          scale_newdata=scale_newdata, center_newdata=center_newdata, realign=realign,
                          scale_reduced=scale_reduced, min_module_genes=min_module_genes)

    return(.se_to_me(newdata, reduced=as.matrix(eig), loadings=loadings(x), assignments=assignments(x), center_X=x@center, scale_X=x@scale))
})

#' @rdname calcEigengenes
#' @export
setMethod("predict", c("ModularExperiment"), function(object, newdata, ...) {
    return(calcEigengenes(object, newdata, ...))
})
