#' FactorisedExperiment: A container for the results of factor analysis
#'
#' @description
#' A container inheriting from \link[ReducedExperiment]{ReducedExperiment}, that
#' contains one or more data matrices, to which factor analysis has been applied
#' to identify a reduced set of features.
#'
#' @param reduced A data matrix, produced by factor analysis, with rows
#' representing samples and columns representing factors.
#'
#' @param loadings A data matrix, produced by factor analysis, with rows
#' representing features and columns representing factors.
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
#' @param stability A vector containing some measure of stability or variance
#' explained for each factor. If factor analysis was performed using
#' \link[ReducedExperiment]{estimate_factors} and `use_stability = True`, this slot will indicate
#' the stability of the factors across multiple runs of ICA.
#'
#' @import SummarizedExperiment
#'
#' @export
FactorisedExperiment <- function(
        loadings = new("matrix"),
        stability = NULL,
        ...)
{
    re <- ReducedExperiment(...)
    return(.FactorisedExperiment(re, loadings=loadings, stability=stability))
}

S4Vectors::setValidity2("FactorisedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Check feature names/numbers
    if (obj_dims[1] != dim(loadings(object))[1])
        msg <- c(msg, "Loadings have invalid row dimensions")

    if (!identical(featureNames(object), rownames(loadings(object))))
        msg <- c(msg, "Loadings have incorrect column names (feature labels)")

    # Factors
    if (dim(loadings(object))[2] != dim(reduced(object))[2])
        msg <- c(msg, "Reduced data and loadings have incompatible column dimensions")

    if (!identical(colnames(loadings(object)), colnames(reduced(object))))
        msg <- c(msg, "Reduced data and loadings have incompatible column names (factor names)")

    # Stability - check names/length matches
    if (!is.null(stability(object)) & length(stability(object)) > 0) {

        if (length(stability(object)) != nComponents(object))
            msg <- c(msg, "Number of components do not match with component stability")

        # If stability vector has names, check they are correct
        if (!is.null(names(stability(object)))) {
            if (!identical(names(stability(object)), componentNames(object)))
                msg <- c(msg, "Component names do not match with component stability")
        }
    }

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("loadings", "FactorisedExperiment", function(object, scale_loadings=FALSE, center_loadings=FALSE, abs_loadings=FALSE) {
    l <- scale(object@loadings, scale=scale_loadings, center=center_loadings)
    if (abs_loadings) l <- abs(l)
    return(l)
})
setReplaceMethod("loadings", "FactorisedExperiment", function(object, value) {
    object@loadings <- value
    validObject(object)
    return(object)
})

setReplaceMethod("names", "FactorisedExperiment", function(x, value) {
    object <- x

    rownames(object@loadings) <- value
    object <- callNextMethod(object, value)
    validObject(object)
    return(object)
})
setReplaceMethod("featureNames", "FactorisedExperiment", function(object, value) {
    names(object) <- value
    return(object)
})
setReplaceMethod("rownames", "FactorisedExperiment", function(x, value) {
    object <- x

    names(object) <- value
    return(object)
})

setMethod("stability", "FactorisedExperiment", function(object) {return(object@stability)})

setReplaceMethod("stability", "FactorisedExperiment", function(object, value) {
    object@stability <- value
    validObject(object)
    return(object)
})

setReplaceMethod("componentNames", "FactorisedExperiment", function(object, value) {
    colnames(object@loadings) <- value
    if (!is.null(object@stability)) names(object@stability) <- value
    object <- callNextMethod(object, value)
    validObject(object)
    return(object)
})

setMethod("[", c("FactorisedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    object <- x

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(object), ",ANY,ANY-method'")

    lod <- object@loadings
    stab <- object@stability

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(object), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(object), fmt
            )
        }
        i <- as.vector(i)
        lod <- lod[i,,drop=FALSE]
    }

    if (!missing(k)) {
        if (is.character(k)) {
            fmt <- paste0("<", class(object), ">[k,] index out of bounds: %s")
            k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                k, componentNames(object), fmt
            )
        }

        k <- as.vector(k)
        lod <- lod[,k,drop=FALSE]
        stab <- stab[k, drop=FALSE]
    }

    out <- callNextMethod(object, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, loadings=lod, stability=stab, check=FALSE)
})

# Same features, different samples
setMethod("cbind", "FactorisedExperiment", function(..., deparse.level=1) {

    args <- list(...)

    loadings_stability_equal <- sapply(args, function(re) {
        return(identical(re@loadings, args[[1]]@loadings) & identical(re@stability, args[[1]]@stability))
    })

    if (!all(loadings_stability_equal))
        stop("Row bind expects loadings and stability slots are equal. Set check_duplicate_slots to FALSE to ignore these slots.")

    args[["deparse.level"]] <- deparse.level

    return(do.call(callNextMethod, args))
})

#' Project new data using pre-defined factors
#'
#' @description
#' Uses a projection approach to calculate factors in new data. Functions in a
#' similar fashion to the `predict` method of \link[stats]{prcomp}. The transposed
#' `newdata` are multiplied by the original loadings matrix.
#'
#' @param object A \link[ReducedExperiment]{FactorisedExperiment} object. The `loadings` slot of
#' this class will be used for projection. Additionally, by default, the `scale`
#' and `center` slots are used to apply the original transformation to the
#' new data.
#'
#' @param newdata New data for projection. Must be a `data.frame` or `matrix`
#' with features as rows and samples as columns, or a
#' \link[SummarizedExperiment]{SummarizedExperiment} object. Assumes that the
#' rows of `newdata` match those of the \link[ReducedExperiment]{FactorisedExperiment}
#' object.
#'
#' @param scale_reduced Whether or not the reduced data should be scaled
#' after calculation.
#'
#' @param scale_newdata Controls whether the `newdata` are scaled. If NULL,
#' performs scaling based on the \link[ReducedExperiment]{FactorisedExperiment}
#' object's `scale` slot. The value of this argument will be passed to the
#' `scale` argument of \link[base]{scale}.
#'
#' @param center_newdata Controls whether the `newdata` are centered If NULL,
#' performs centering based on the \link[ReducedExperiment]{FactorisedExperiment}
#' object's `center` slot. The value of this argument will be passed to the
#' `center` argument of \link[base]{scale}.
#'
#' @param assay_name If a \link[SummarizedExperiment]{SummarizedExperiment}
#' object is passed as new data, this argument indicates which assay should be
#' used for projection.
#'
#' @returns Calculates a matrix with samples as rows and factors as columns. If
#' `newdata` was a `matrix` or `data.frame`, this will be returned as a matrix.
#' If a \link[SummarizedExperiment]{SummarizedExperiment} object was passed
#' instead, then a If a \link[ReducedExperiment]{FactorisedExperiment}
#' object will be created containing this matrix in its `reduced` slot.
#'
#' @seealso [ReducedExperiment::calcEigengenes()]
#'
#' @rdname projectData
#' @export
setMethod("projectData", c("FactorisedExperiment", "matrix"), function(object, newdata, scale_reduced=TRUE, scale_newdata=NULL, center_newdata=NULL) {

    if (!identical(rownames(object), rownames(newdata)))
        stop("Rownames of x do not match those of newdata")

    # apply known vectors for scaling and centering (returned as attributes by `scale`)
    if (is.null(scale_newdata)) scale_newdata <- object@scale
    if (is.null(center_newdata)) center_newdata <- object@center

    newdata <- t(scale(t(newdata), scale=scale_newdata, center=center_newdata))
    red <- .project_ica(newdata, loadings(object))

    if (scale_reduced) red <- scale(red)

    return(red)
})

#' @rdname projectData
#' @export
setMethod("projectData", c("FactorisedExperiment", "data.frame"), function(object, newdata, scale_reduced=TRUE, scale_newdata=NULL, center_newdata=NULL) {
    return(projectData(object, as.matrix(newdata), scale_reduced=scale_reduced, scale_newdata=scale_newdata, center_newdata=center_newdata))
})

#' @rdname projectData
#' @export
setMethod("projectData", c("FactorisedExperiment", "SummarizedExperiment"), function(object, newdata, scale_reduced=TRUE, scale_newdata=NULL, center_newdata=NULL, assay_name="normal") {

    projected_data <- projectData(object, assay(newdata, assay_name), scale_reduced=scale_reduced, scale_newdata=scale_newdata, center_newdata=center_newdata)

    return(.se_to_fe(newdata, reduced=projected_data, loadings=loadings(object), stability=stability(object), center_X=object@center, scale_X=object@scale))
})

#' @rdname projectData
#' @export
setMethod("predict", c("FactorisedExperiment"), function(object, newdata, ...) {
    return(projectData(object, newdata, ...))
})

setMethod("getAlignedFeatures", c("FactorisedExperiment"), function(object, z_cutoff=NULL, n_features=NULL,
                                                                    feature_id_col="rownames", format="list",
                                                                    scale_loadings=TRUE, center_loadings=FALSE) {

    S <- loadings(object, scale_loadings=scale_loadings, center_loadings=center_loadings)

    if (feature_id_col != "rownames") rownames(S) <- rowData(object)[[feature_id_col]]

    if (is.null(z_cutoff) & is.null(n_features)) n_features <- nrow(S)
    if (is.null(n_features)) n_features <- 0
    if (n_features > nrow(S)) n_features <- nrow(S)

    factor_features <- data.frame()
    for (f in componentNames(object)) {

        which_features <- c()
        if (!is.null(z_cutoff)) which_features <- which(abs(S[,f]) > z_cutoff)

        if (length(which_features) < n_features) {

            which_features <- order(abs(S[,f]), decreasing = TRUE)[1:n_features]
        }

        if (length(which_features) > 0) {
            factor_features <- rbind(factor_features, data.frame(
                component = f,
                feature = rownames(S)[which_features],
                value = S[,f][which_features],
                loadings_scaled = scale_loadings,
                loadings_centered = center_loadings
            ))
        }
    }

    factor_features$z_cutoff = z_cutoff
    factor_features$n_features = n_features

    if (is.function(format)) {
        return(format(factor_features))
    } else if (format == "data.frame") {
        return(factor_features)
    } else if (format == "list") {

        factor_list <- list()

        for (f in unique(factor_features$component))
            factor_list[[f]] <- factor_features$feature[which(factor_features$component == f)]

        return(factor_list)
    }
})

setMethod("runEnrich", c("FactorisedExperiment"),
    function(object, method="overrepresentation", feature_id_col="rownames",
             scale_loadings=TRUE, abs_loadings=FALSE, z_cutoff=3, n_features=20,
             as_dataframe=FALSE, ...)
{
    if (method == "gsea") {

        z_cutoff <- NULL
        n_features <- NULL

        S <- loadings(object, scale_loadings=scale_loadings, abs_loadings=abs_loadings)
        if (feature_id_col != "rownames") rownames(S) <- rowData(object)[[feature_id_col]]

        enrich_res <- reduced_gsea(S, ...)

    } else if (method == "overrepresentation") {

        factor_features <- getAlignedFeatures(object, feature_id_col=feature_id_col,
                                              scale_loadings=scale_loadings,
                                              z_cutoff=z_cutoff,
                                              n_features=n_features)

        enrich_res <- reduced_oa(factor_features, ...)

    } else {
        stop("Enrichment method not recognised")
    }

    for (comp in names(enrich_res)) {
        if (!is.null(enrich_res[[comp]]@result)) {
            if (nrow(enrich_res[[comp]]@result) >= 1) {
                enrich_res[[comp]]@result$z_cutoff <- z_cutoff
                enrich_res[[comp]]@result$n_features <- n_features
                enrich_res[[comp]]@result$loadings_scaled <- scale_loadings
                enrich_res[[comp]]@result$loadings_centered <- FALSE
                enrich_res[[comp]]@result$abs_loadings <- abs_loadings
            }
        }
    }

    if (as_dataframe)  {
        enrich_res <- lapply(enrich_res, function(object) {object@result})
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})
