#' FactorisedExperiment
#'
#' A class to represent the results of factor analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
FactorisedExperiment <- function(
        loadings = new("matrix"),
        stability = new("numeric"),
        scale = FALSE,
        center = FALSE,
        ...)
{
    re <- ReducedExperiment(...)
    return(.FactorisedExperiment(re, loadings=loadings, stability=stability, scale=scale, center=center))
}

S4Vectors::setValidity2("FactorisedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Check feature names
    if (obj_dims[1] != dim(loadings(object))[1])
        msg <- c(msg, "Loadings have invalid row dimensions")
    if (!all.equal(featureNames(object), rownames(loadings(object))))
        msg <- c(msg, "Loadings have incorrect column names (feature labels)")

    # Factors
    if (dim(loadings(object))[2] != dim(reduced(object))[2])
        msg <- c(msg, "Reduced data and loadings have incompatible column dimensions")

    if (!all.equal(colnames(loadings(object)), colnames(reduced(object))))
        msg <- c(msg, "Reduced data and loadings have incompatible column names (factor names)")

    # Center / Scale
    if (!is.logical(object@scale)) {
        if (!all.equal(rownames(object), names(object@scale)))
            msg <- c(msg, "Scaling vector has invalid names")
    }
    if (!is.logical(object@scale)) {
        if (!all.equal(rownames(object), names(object@center)))
            msg <- c(msg, "Centering vector has invalid names")
    }

    # Stability
    if (!is.null(stability(object))) {
        if (length(stability(object)) != nComponents(object))
            msg <- c(msg, "Number of components do not match with component stability")

        if (!all.equal(names(stability(object)), componentNames(object)))
            msg <- c(msg, "Component names do not match with component stability")
    }

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("loadings", "FactorisedExperiment", function(x, scale_X=FALSE, center_X=FALSE) {
    return(scale(x@loadings, scale=scale_X, center=center_X))
})

setReplaceMethod("loadings", "FactorisedExperiment", function(x, value) {
    x@loadings <- value
    validObject(x)
    return(x)
})

setMethod("stability", "FactorisedExperiment", function(x) {return(x@stability)})

setReplaceMethod("stability", "FactorisedExperiment", function(x, value) {
    x@stability <- value
    validObject(x)
    return(x)
})

setReplaceMethod("componentNames", "FactorisedExperiment", function(x, value) {
    colnames(x@loadings) <- value
    x <- callNextMethod(x, value)
    validObject(x)
    return(x)
})

setMethod("[", c("FactorisedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    lod <- loadings(x)
    stab <- x@stability
    center <- x@center
    scale <- x@scale

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        lod <- lod[i,,drop=FALSE]
        if (!is.logical(center)) center <- center[i,drop=FALSE]
        if (!is.logical(scale)) scale <- scale[i,drop=FALSE]
    }

    if (!missing(k)) {
        if (is.character(k)) {
            fmt <- paste0("<", class(x), ">[k,] index out of bounds: %s")
            k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                k, componentNames(x), fmt
            )
        }

        k <- as.vector(k)
        lod <- lod[,k,drop=FALSE]
        stab <- stab[k, drop=FALSE]
    }

    out <- callNextMethod(x, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, loadings=lod, stability=stab, center=center, scale=scale, check=FALSE)
})

#' Project data
setMethod("projectData", c("FactorisedExperiment", "matrix"), function(x, newdata) {

    newdata <- t(scale(t(newdata), scale=x@scale, center=x@center))

    return(.project_ica(newdata, loadings(x)))
})

setMethod("projectData", c("FactorisedExperiment", "data.frame"), function(x, newdata) {
    return(projectData(x, as.matrix(newdata)))
})

setMethod("projectData", c("FactorisedExperiment", "SummarizedExperiment"), function(x, newdata, assay_name="normal") {

    projected_data <- projectData(x, assay(newdata, assay_name))

    return(.se_to_fe(newdata, reduced=projected_data, loadings=loadings(x), stability=stability(x), center=x@center, scale=x@scale))
})

setMethod("predict", c("FactorisedExperiment"), function(object, newdata, ...) {
    return(projectData(object, newdata, ...))
})

setMethod("getAlignedFeatures", c("FactorisedExperiment"), function(x, z_cutoff=NULL, n_features=NULL,
                                                                    feature_id_col="rownames", format="list",
                                                                    scale_loadings=TRUE) {
    S <- loadings(x, scale=scale_loadings)
    if (feature_id_col != "rownames") rownames(S) <- rowData(x)[[feature_id_col]]

    if (is.null(z_cutoff) & is.null(n_features)) n_features <- nrow(S)
    if (n_features > nrow(S)) n_features <- nrow(S)

    factor_features <- data.frame()
    for (f in componentNames(x)) {

        which_features <- c()
        if (!is.null(z_cutoff)) which_features <- which(abs(S[,f]) > z_cutoff)

        if (length(which_features) < n_features) {

            which_features <- order(abs(S[,f]), decreasing = TRUE)[1:n_features]
        }

        factor_features <- rbind(factor_features, data.frame(
            component = f,
            feature = rownames(S)[which_features],
            value = S[,f][which_features],
            loadings_scaled = scale_loadings,
            loadings_centered = TRUE
        ))
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
    function(x, method="overrepresentation", feature_id_col="rownames",
             scale_loadings=TRUE, z_cutoff=3, n_features=20, as_dataframe=FALSE, ...)
{
    if (method == "gsea") {

        z_cutoff <- NULL
        n_features <- NULL

        S <- loadings(x, scale=scale_loadings)
        if (feature_id_col != "rownames") rownames(S) <- rowData(x)[[feature_id_col]]

        enrich_res <- reduced_gsea(S, ...)

    } else if (method == "overrepresentation") {

        factor_features <- getAlignedFeatures(x, feature_id_col=feature_id_col,
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
            }
        }
    }

    if (as_dataframe)  {
        enrich_res <- lapply(enrich_res, function(x) {x@result})
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})
