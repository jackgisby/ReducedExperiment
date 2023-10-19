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

setMethod("loadings", "ReducedExperiment", function(x, withDimnames=TRUE) {
    out <- x@loadings
    if (withDimnames) {
        rownames(out) <- rownames(x)
        colnames(out) <- componentNames(x)
    }
    return(out)
})

setReplaceMethod("loadings", "ReducedExperiment", function(x, value) {
    x@loadings <- value
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
          function(x, i, j, k, ..., drop=FALSE) {

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    lod <- loadings(x, withDimnames=FALSE)
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
setMethod("project_data", c("FactorisedExperiment", "matrix"), function(x, newdata) {
    return(.project_ica(newdata, loadings(x)))
})

setMethod("project_data", c("FactorisedExperiment", "SummarizedExperiment"), function(x, newdata) {
    return(.project_ica(assay(newdata, "normal"), loadings(x)))
})

# setMethod("predict", signature="FactorisedExperiment",  function(object, newdata) {
#     # TODO: Implement function
#     stop("Not implemented yet.")
#     return(newdata)
# })
