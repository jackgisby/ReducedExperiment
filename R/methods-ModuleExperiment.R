#' ModularExperiment
#'
#' A class to represent the results of module analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
ModularExperiment <- function(assignments=character(), ...)
{
    re <- ReducedExperiment(...)
    return(.ModularExperiment(re, assignments=assignments))
}

S4Vectors::setValidity2("ModularExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Features
    if (obj_dims[3] != length(assignments(object)))
        msg <- c(msg, "Assignments have invalid length")

    if (!all.equal(assignments(object), rownames(object)))
        msg <- c(msg, "Assignments have incompatible names (rownames)")

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

setMethod("moduleNames", "ModularExperiment", function(x, as_list=FALSE) {
    return(assignments(x, as_list))
})

setReplaceMethod("moduleNames", "ModularExperiment", function(x, value) {
    assignments(x) <- value
    return(x)
})

setReplaceMethod("componentNames", "ModularExperiment", function(x, value) {
    names(x@assignments) <- value
    x <- callNextMethod(x, value)
    validObject(x)
    return(x)
})

setMethod("[", c("ModularExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    assignments <- x@assignments

    if (!missing(k)) {
        if (is.character(k)) {
            fmt <- paste0("<", class(x), ">[k,] index out of bounds: %s")
            k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
              k, componentnames(x), fmt
            )
      }

        k <- as.vector(k)
        assignments <- assignments[k, drop=FALSE]
    }

    out <- callNextMethod(x, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, assignments=assignments, check=FALSE)
})

setMethod("runEnrich", c("ModularExperiment"),
          function(x, method="overrepresentation", feature_id_col="rownames", as_dataframe=FALSE, ...)
{
    if (method == "overrepresentation") {

        modules <- assignments(x)

        enrich_res <- reduced_oa(modules, ...)

    } else {
        stop("Enrichment method not recognised")
    }

    # for (comp in names(enrich_res)) {
    #     if (!is.null(enrich_res[[comp]]@result)) {
    #         if (nrow(enrich_res[[comp]]@result) >= 1) {
    #           enrich_res[[comp]]@result$setting <- TRUE
    #         }
    #     }
    # }

    if (as_dataframe) enrich_res <- do.call("rbind", enrich_res)

    return(enrich_res)
})
