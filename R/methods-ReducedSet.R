#' ReducedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
ReducedExperiment <- function(
        reduced = new("matrix"),
        ...)
{
    se <- SummarizedExperiment::SummarizedExperiment(...)
    return(.ReducedExperiment(se, reduced=reduced))
}

S4Vectors::setValidity2("ReducedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Samples
    if (obj_dims[2] != dim(reduced(object))[1])
        msg <- c(msg, "Reduced data have invalid row dimensions")

    rnames <- if (is.null(rownames(reduced(object)))) character(0) else rownames(reduced(object))
    if (!all.equal(sampleNames(object), rnames))
        msg <- c(msg, "Reduced data have invalid row names (sample labels)")

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("reduced", "ReducedExperiment", function(x, withDimnames=TRUE) {
    out <- x@reduced
    if (withDimnames) {
        rownames(out) <- colnames(x)
        colnames(out) <- componentNames(x)
    }
    return(out)
})

setReplaceMethod("reduced", "ReducedExperiment", function(x, value) {
    x@reduced <- value
    validObject(x)
    return(x)
})

setMethod("componentNames", "ReducedExperiment",
          function(x) {return(colnames(x@reduced))})

setReplaceMethod("componentNames", "ReducedExperiment", function(x, value) {
    colnames(x@reduced) <- value
    validObject(x)
    return(x)
})

# TODO: setters
setMethod("featureNames", "ReducedExperiment", function(x) {return(names(x))})
setMethod("sampleNames", "ReducedExperiment", function(x) {return(rownames(colData(x)))})

setMethod("show", "ReducedExperiment" ,
          function(object) {
              callNextMethod()
              cat(nComponents(object), " components\n")
          })

setMethod("[", c("ReducedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE) {

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    red <- reduced(x, withDimnames=FALSE)

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                j, colnames(x), fmt
            )
        }
        j <- as.vector(j)
        red <- red[j,,drop=FALSE]
    }


    if (!missing(k)) {
        if (is.character(k)) {
            fmt <- paste0("<", class(x), ">[k,] index out of bounds: %s")
            k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                k, componentnames(x), fmt
            )
        }
        k <- as.vector(k)
        red <- red[,k,drop=FALSE]
    }

    out <- callNextMethod(x, i, j, ...)
    BiocGenerics:::replaceSlots(out, reduced=red, check=FALSE)
})

setMethod("dim", "ReducedExperiment", function(x) {
    out <- c(callNextMethod(x), ncol(reduced(x)))
    names(out) <- c("Features", "Samples", "Components")
    return(out)
})

setMethod("nComponents", "ReducedExperiment", function(x) {dim(x)[3]})
setMethod("nSamples", "ReducedExperiment", function(x) {dim(x)[2]})
setMethod("nFeatures", "ReducedExperiment", function(x) {dim(x)[1]})
