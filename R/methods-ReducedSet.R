setMethod(
    "initialize",
    signature = "ReducedSet",
    definition = function(.Object,
                          assayData,
                          phenoData,
                          featureData,
                          reducedData,
                          L,
                          exprs = new("matrix"),
                          ...)
    {
        #TODO: store X scaling/centering attr (in assayData or own slot)
        #TODO: store original S / M but provide options to scale/center

        # New slots
        .Object@reducedData <- reducedData
        .Object@L <- L

        # For consistency with ExpressionSet
        if (missing(assayData)) {
            if (missing(phenoData))
                phenoData <-
                    Biobase::annotatedDataFrameFrom(exprs, byrow = FALSE)
            if (missing(featureData))
                featureData <-
                    Biobase::annotatedDataFrameFrom(exprs, byrow = TRUE)
            .Object <- callNextMethod(
                .Object,
                phenoData = phenoData,
                featureData = featureData,
                exprs = exprs,
                ...
            )
        } else {
            if (missing(phenoData))
                phenoData <-
                    Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE)
            if (missing(featureData))
                featureData <-
                    Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE)
            .Object <- callNextMethod(
                .Object,
                assayData = assayData,
                phenoData = phenoData,
                featureData = featureData,
                ...
            )
        }

        Biobase:::.harmonizeDimnames(.Object)
    }
)

setValidity("ReducedSet", function(object) {
    msg <-
        Biobase::validMsg(NULL, Biobase:::isValidVersion(object, "ReducedSet"))
    msg <-
        validMsg(msg, Biobase::assayDataValidMembers(assayData(object), "exprs"))

    obj_dims <- dim(object)
    # Features
    if (obj_dims[1] != dim(loadings(object))[1])
        msg <-
        Biobase::validMsg(msg, "Loadings have invalid row dimensions")
    if (!all.equal(featureNames(object), rownames(loadings(object))))
        msg <-
        Biobase::validMsg(msg, "Loadings have incorrect column names (feature labels)")

    # Samples
    if (obj_dims[2] != dim(reduced(object))[1])
        msg <-
        Biobase::validMsg(msg, "Reduced data have invalid row dimensions")
    if (!all.equal(sampleNames(object), rownames(reduced(object))))
        msg <-
        Biobase::validMsg(msg, "Reduced data have invalid row names (sample labels)")

    # Metagenes
    if (dim(loadings(object))[2] != dim(reduced(object))[2])
        msg <-
        Biobase::validMsg(msg,
                          "Reduced data and loadings have incompatible column dimensions")

    if (!all.equal(colnames(loadings(object)), colnames(reduced(object))))
        msg <-
        Biobase::validMsg(
            msg,
            "Reduced data and loadings have incompatible column names (module/factor names)"
        )

    return(if (is.null(msg))
        TRUE
        else
            msg)
})

setAs("ReducedSet", "data.frame",
      function (from)
          data.frame(reduced(from), pData(from)))

as.data.frame.ReducedSet <- Biobase:::as.data.frame.ExpressionSet

setMethod("reducedData", "ReducedSet", function(object)
    object@reducedData)

setReplaceMethod("reducedData",
                 signature = signature(object = "ReducedSet",
                                       value = "AssayData"),
                 function(object, value) {
                     object@reducedData <- value
                     object
                 })

setMethod("exprs", signature(object = "ReducedSet"),
          function(object)
              assayDataElement(object, "exprs"))

setReplaceMethod("exprs", signature(object = "ReducedSet", value = "matrix"),
                 function(object, value)
                     assayDataElementReplace(object, "exprs", value))

setMethod("reduced", signature(object = "ReducedSet"),
          function(object)
              object@reducedData)

setReplaceMethod("reduced", signature(object = "ReducedSet", value = "matrix"),
                 function(object, value) {
                     object@reducedData <- value
                     object
                 })

setMethod("loadings", signature(object = "ReducedSet"), function(object)
    object@L)

setReplaceMethod("loadings", signature(object = "ReducedSet", value = "matrix"),
                 function(object, value) {
                     object@L <- value
                     return(object)
                 })

setMethod("write.exprs",
          signature(x = "ReducedSet"),
          function(x,
                   file = "tmp.txt",
                   quote = FALSE,
                   sep = "\t",
                   col.names = NA,
                   ...) {
              write.table(
                  exprs(x),
                  file = file,
                  quote = quote,
                  sep = sep,
                  col.names = col.names,
                  ...
              )
          })

setMethod("write.reduced",
          signature(x = "ReducedSet"),
          function(x,
                   file = "tmp.txt",
                   quote = FALSE,
                   sep = "\t",
                   col.names = NA,
                   ...) {
              write.table(
                  reduced(x),
                  file = file,
                  quote = quote,
                  sep = sep,
                  col.names = col.names,
                  ...
              )
          })

setMethod("show", "ReducedSet" ,
          function(object) {
              cat(nDim(object), " components\n")
              callNextMethod()
          })

setMethod("[",
          "ReducedSet",
          function(x, i, j, k, ..., drop = FALSE)
          {
              if (!missing(k)) {

              }
          })
