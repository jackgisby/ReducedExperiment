setMethod(
    "initialize",
    signature="ReducedSet",
    definition=function(
        .Object,
        assayData=Biobase::assayDataNew(exprs=exprs),
        reducedData=new("matrix"),
        phenoData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
        featureData=Biobase::annotatedDataFrameFrom(assayData, byrow=TRUE),
        exprs=new("matrix"),
        ...)
{
    if (!missing(exprs) & !missing(assayData)) {
        stop("Only specify one of `exprs` and `assayData`")
    }

    #TODO: store X scaling/centering attr (in assayData or own slot)

    # Reduced data common to children
    .Object@reducedData <- reducedData

    callNextMethod(.Object, assayData=assayData, phenoData=phenoData,
                   featureData=featureData, ...)
})

setValidity("ReducedSet", function(object) {
    msg <- Biobase::validMsg(
        NULL,
        Biobase:::isValidVersion(object, "ReducedSet")
    )

    obj_dims <- dim(object)

    # Samples
    if (obj_dims[2] != dim(reduced(object))[1])
        msg <- Biobase::validMsg(
            msg,
            "Reduced data have invalid row dimensions"
        )

    rnames <- if (is.null(rownames(reduced(object)))) character(0) else rownames(reduced(object))
    if (!all.equal(sampleNames(object), rnames))
        msg <- Biobase::validMsg(
            msg,
            "Reduced data have invalid row names (sample labels)"
        )

    return(if (is.null(msg)) TRUE else msg)
})

setAs("ReducedSet", "data.frame",
      function (from) {data.frame(reduced(from), pData(from))})

as.data.frame.ReducedSet <- Biobase:::as.data.frame.ExpressionSet

setMethod("reduced", signature(object="ReducedSet"),
          function(object) {return(object@reducedData)})

setReplaceMethod("reduced",
                 signature=signature(object="ReducedSet", value="matrix"),
                 function(object, value) {
                     object@reducedData <- value
                     return(object)
                 })

setMethod("exprs", signature(object="ReducedSet"),
          function(object) Biobase::assayDataElement(object,"exprs"))

setReplaceMethod("exprs", signature(object="ReducedSet",value="matrix"),
                 function(object, value) Biobase::assayDataElementReplace(object, "exprs", value))

setMethod("write.reduced",
          signature(x="ReducedSet"),
          function(x, file="tmp.txt", quote=FALSE, sep="\t", col.names=NA,
                   ...) {
              write.table(reduced(x), file=file, quote=quote, sep=sep,
                          col.names=col.names, ...)
          })

setMethod("show", "ReducedSet" ,
          function(object) {
              cat(nComponents(object), " components\n")
              return(callNextMethod())
          })

.whichToKeep <- function(idx, xnames) {
    if (is.character(idx)) {
        x_to_keep <- which(xnames %in% idx)
    } else if (is.numeric(idx)) {
        x_to_keep <- idx
    } else if (is.logical(idx)) {
        x_to_keep <- which(idx)
    } else {
        stop("Provide a character, numeric or logical vector to slice.")
    }

    return(x_to_keep)
}

setMethod("[", "ReducedSet",
          function(x, i, j, k, ..., drop=FALSE)
{
    if (missing(i) & missing(j) & missing(k)) {
        stop("Specify a dimension for slicing")
    }

    if (!missing(i) | !missing(j)) {
        x <- callNextMethod(x, i, j, ..., drop=drop)

        if (!missing(j)) {  # Samples
            reduced(x) <- reduced(x)[.whichToKeep(j, rownames(reduced(x))), ]
        }
    }

    if (!missing(k)) {  # Components
        reduced(x) <- reduced(x)[, .whichToKeep(k, colnames(reduced(x)))]
    }

    return(x)
})

setMethod("dim", "ReducedSet", function(x)
    {c(callNextMethod(x), "Components"=ncol(reduced(x)))})

setMethod("dims", "ReducedSet", function(x)
    {rbind(callNextMethod(x), matrix(ncol(reduced(x)),  dimnames=list("Components")))})

setMethod("nComponents", "ReducedSet", function(object) {dim(object)[3]})
setMethod("nSamples", "ReducedSet", function(object) {dim(object)[2]})
setMethod("nFeatures", "ReducedSet", function(object) {dim(object)[1]})
