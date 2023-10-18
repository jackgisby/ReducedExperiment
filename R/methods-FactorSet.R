setMethod(
    "initialize",
    signature="FactorSet",
    definition=function(.Object, S, varexp=new("numeric"),
                        assayData=Biobase::assayDataNew(exprs=exprs),
                        reducedData=new("matrix"),
                        phenoData=Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE),
                        featureData=Biobase::annotatedDataFrameFrom(assayData, byrow=TRUE),
                        exprs=new("matrix"), ...) {

        if (!missing(exprs) & !missing(assayData)) {
            stop("Only specify one of `exprs` and `assayData`")
        }

        .Object@S <- S
        .Object@varexp <- varexp

        callNextMethod(.Object, assayData=assayData, reducedData=reducedData, phenoData=phenoData, featureData=featureData, ...)
    }
)    #TODO: store original S / M but provide options to scale/center

setValidity("FactorSet", function(object) {
    msg <- Biobase::validMsg(
        NULL,
        Biobase:::isValidVersion(object, "FactorSet")
    )

    obj_dims <- dim(object)

    # Features
    if (obj_dims[1] != dim(loadings(object))[1])
        msg <- Biobase::validMsg(msg, "Loadings have invalid row dimensions")
    if (!all.equal(featureNames(object), rownames(loadings(object))))
        msg <- Biobase::validMsg(
            msg,
            "Loadings have incorrect column names (feature labels)"
        )

    # Factors
    if (dim(loadings(object))[2] != dim(reduced(object))[2])
        msg <- Biobase::validMsg(
            msg,
            "Reduced data and loadings have incompatible column dimensions"
        )

    if (!all.equal(colnames(loadings(object)), colnames(reduced(object))))
        msg <- Biobase::validMsg(
            msg,
            "Reduced data and loadings have incompatible column names (factor
            names)"
        )

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("loadings", signature(object="FactorSet"),
          function(object) {return(object@S)})

setReplaceMethod("loadings",
                 signature=signature(object="FactorSet", value="matrix"),
                 function(object, value) {
                     object@S <- value
                     return(object)
                 })

setMethod("[", "FactorSet",
          function(x, i, j, k, ..., drop=FALSE)
{
    x <- callNextMethod(x, i, j, k, ..., drop=drop)

    if (!missing(i)) {  # Features
      loadings(x) <- loadings(x)[.whichToKeep(i, rownames(loadings(x))), ]
    }

    if (!missing(k)) {  # Components
      loadings(x) <- loadings(x)[, .whichToKeep(k, colnames(reduced(x)))]
    }

    return(x)
})
