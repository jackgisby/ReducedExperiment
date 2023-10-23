#' FactorisedExperiment
#'
#' A class to represent the results of factor analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
FactorisedExperiment <- function(
        loadings = new("matrix"),
        varexp = new("numeric"),
        ...)
{
    re <- ReducedExperiment(...)
    return(.FactorisedExperiment(re, loadings=loadings, varexp=varexp))
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

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                i, rownames(x), fmt
            )
        }
        i <- as.vector(i)
        lod <- lod[i,,drop=FALSE]
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
    BiocGenerics:::replaceSlots(out, loadings=lod, varexp=vafs, check=FALSE)
})

#' Project data
setMethod("project_data", c("FactorisedExperiment", "matrix"), function(x, newdata, new_data_is_transformed=FALSE) {

    if (!new_data_is_transformed) {
        newdata <- .scale_center_based_on_attr(assay(x, "transformed"), newdata)
    }

    return(.project_ica(newdata, loadings(x)))
})

setMethod("project_data", c("FactorisedExperiment", "SummarizedExperiment"), function(x, newdata, new_data_is_transformed=FALSE) {

    if (!new_data_is_transformed) {
        assay(newdata, "transformed") <- .scale_center_based_on_attr(assay(x, "transformed"), assay(newdata, "normal"))
    }

    projected_data <- project_data(x, assay(newdata, "transformed"), new_data_is_transformed=TRUE)

    return(.se_to_fe(newdata, reduced=projected_data, loadings=loadings(x), varexp=varexp(x)))
})

.scale_center_based_on_attr <- function(x, newdata) {
    if ("scaled:center" %in% names(attributes(x))) {
        center <- attr(x, "scaled:center")
    } else {
        center <- FALSE
    }
    if ("scaled:scale" %in% names(attributes(x))) {
        scale <- attr(x, "scaled:scale")
    } else {
        scale <- FALSE
    }

    return(t(scale(t(newdata), scale=scale, center=center)))
}

setMethod("predict", c("FactorisedExperiment"), function(object, newdata, new_data_is_transformed=FALSE) {
    return(project_data(object, newdata, new_data_is_transformed))
})

setMethod("run_enrichment", c("FactorisedExperiment"),
    function(object, ensembl_id_col="ensembl_id", entrezgene_id=NULL, scale_loadings=TRUE) {

        factor_genes <- list()
        for (f in componentNames(object)) {

            factor_genes[[f]] <-
            enrich_res <- rbind(enrich_res, comp_enrich)
        }
})
