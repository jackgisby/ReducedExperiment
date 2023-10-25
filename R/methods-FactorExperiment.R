#' FactorisedExperiment
#'
#' A class to represent the results of factor analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
FactorisedExperiment <- function(
        loadings = new("matrix"),
        varexp = new("numeric"),
        center=FALSE,
        ...)
{
    re <- ReducedExperiment(...)
    return(.FactorisedExperiment(re, loadings=loadings, varexp=varexp, center=center))
}

S4Vectors::setValidity2("FactorisedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Features
    if (obj_dims[1] != dim(loadings(object))[1])
        msg <- c(msg, "Loadings have invalid row dimensions")
    if (!all.equal(featureNames(object), rownames(loadings(object))))
        msg <- c(msg, "Loadings have incorrect column names (feature labels)")

    # Factors
    if (dim(loadings(object))[2] != dim(reduced(object))[2])
        msg <- c(msg, "Reduced data and loadings have incompatible column dimensions")

    if (!all.equal(colnames(loadings(object)), colnames(reduced(object))))
        msg <- c(msg, "Reduced data and loadings have incompatible column names (factor names)")

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("loadings", "ReducedExperiment", function(x, scale=FALSE, center=FALSE) {
    return(scale(x@loadings, scale=scale, center=center))
})

setReplaceMethod("loadings", "ReducedExperiment", function(x, value) {
    x@loadings <- value
    validObject(x)
    return(x)
})

setMethod("varexp", "ReducedExperiment", function(x) {return(x@varexp)})

setReplaceMethod("varexp", "ReducedExperiment", function(x, value) {
    x@varexp <- value
    validObject(x)
    return(x)
})

setReplaceMethod("componentNames", "ReducedExperiment", function(x, value) {
    colnames(x@reduced) <- value
    colnames(x@loadings) <- value
    validObject(x)
    return(x)
})

setMethod("[", c("FactorisedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    lod <- loadings(x)
    vafs <- x@varexp
    center <- x@center

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        lod <- lod[i,,drop=FALSE]
        center <- center[i,drop=FALSE]
    }

    if (!missing(k)) {
        if (is.character(k)) {
            fmt <- paste0("<", class(x), ">[k,] index out of bounds: %s")
            k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                k, componentnames(x), fmt
            )
        }

        k <- as.vector(k)
        lod <- lod[,k,drop=FALSE]
        vafs <- vafs[k, drop=FALSE]
    }

    out <- callNextMethod(x, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, loadings=lod, varexp=vafs, center=center, check=FALSE)
})

#' Project data
setMethod("project_data", c("FactorisedExperiment", "matrix"), function(x, newdata, new_data_is_centered=FALSE) {

    if (!new_data_is_centered) newdata <- t(scale(t(newdata), scale=FALSE, center=x@center))

    return(.project_ica(newdata, loadings(x)))
})

setMethod("project_data", c("FactorisedExperiment", "SummarizedExperiment"), function(x, newdata, new_data_is_centered=FALSE) {

    if (!new_data_is_centered) {
        assay(newdata, "transformed") <- t(scale(t(assay(newdata, "normal")), scale=FALSE, center=x@center))
    }

    projected_data <- project_data(x, assay(newdata, "transformed"), new_data_is_centered=TRUE)

    return(.se_to_fe(newdata, reduced=projected_data, loadings=loadings(x), varexp=varexp(x), center=x@center))
})

setMethod("predict", c("FactorisedExperiment"), function(object, newdata, new_data_is_transformed=FALSE) {
    return(project_data(object, newdata, new_data_is_transformed))
})

setMethod("runEnrich", c("FactorisedExperiment"),
    function(x, method=c("overrepresentation", "gsea"), feature_id_col="rownames",
             scale_loadings=TRUE, z_cutoff=3, min_features=20, ...)
{
    S <- loadings(x, scale=scale_loadings)
    if (feature_id_col != "rownames") rownames(S) <- rowData(x)[[feature_id_col]]

    if (method == "gsea") {
        enrich_res <- reduced_gsea(S, ...)
    }

    factor_features <- list()
    for (f in componentNames(x)) {

        factor_features[[f]] <- rownames(S)[which(abs(S[,f]) > z_cutoff)]

        if (length(factor_features[[f]]) < min_features) {

            factor_features[[f]] <- rownames(S)[order(abs(S[,f]), decreasing = TRUE)][1:20]
        }
    }

    if (method == "overrepresentation") {
        enrich_res <- reduced_oa(factor_features, ...)
    } else {
        stop(paste0("Enrichment method ", method, " is not a valid option"))
    }

    if (nrow(enrich_res) >= 1) {
        enrich_res$z_cutoff <- z_cutoff
        enrich_res$min_features <- min_features
        enrich_res$loadings_scaled <- scale_loadings

        enrich_res <- subset(enrich_res, select=c("ID", "description", "factor",
            "z_cutoff", "min_features", "loadings_scaled", "gene_ratio",
            "bg_ratio", "pvalue", "method", "adj_pvalue","adj_method", "count",
            "gene_id"))
    }

    return(enrich_res)
})
