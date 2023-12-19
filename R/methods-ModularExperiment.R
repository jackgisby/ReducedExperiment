#' ModularExperiment
#'
#' A class to represent the results of module analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
ModularExperiment <- function(assignments=character(), dendrogram=NULL, threshold=NULL, ...)
{
    re <- ReducedExperiment(...)
    return(.ModularExperiment(re, assignments=assignments, dendrogram=dendrogram, threshold=threshold))
}

S4Vectors::setValidity2("ModularExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Assignments
    if (obj_dims[1] != length(assignments(object)))
        msg <- c(msg, "Assignments have invalid length")

    if (!all(assignments(object) == rownames(object)))
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

setReplaceMethod("featureNames", "ModularExperiment", function(x, value) {
    callNextMethod(x, value)
    x@assignments <- setNames(x@assignments, names(x))
    validObject(x)
    return(x)
})

setReplaceMethod("componentNames", "ModularExperiment", function(x, value) {

    curr_names <- colnames(x@reduced)
    x <- callNextMethod(x, value)
    new_names <- colnames(x@reduced)

    for (i in 1:length(curr_names)) {
        names(x@assignments)[which(names(x@assignments) == curr_names[i])] <- new_names[i]
    }

    validObject(x)
    return(x)
})

setMethod("moduleNames", "ModularExperiment", function(x) {
    return(componentNames(x))
})

setReplaceMethod("moduleNames", "ModularExperiment", function(x, value) {
    componentNames(x) <- value
    return(x)
})

setMethod("nModules", "ModularExperiment", function(x) {dim(x)[3]})

setMethod("dendrogram", "ModularExperiment", function(x) {
    return(x@dendrogram)
})

setReplaceMethod("dendrogram", "ModularExperiment", function(x, value) {
    dendrogram(x) <- value
    return(x)
})

setMethod("[", c("ModularExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    assignments <- x@assignments

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
              i, rownames(x), fmt
            )
      }

        i <- as.vector(i)
        assignments <- assignments[i, drop=FALSE]
    }

    out <- callNextMethod(x, i, j, k, ...)
    BiocGenerics:::replaceSlots(out, assignments=assignments, check=FALSE)
})

setMethod("runEnrich", c("ModularExperiment"),
          function(x, method="overrepresentation", feature_id_col="rownames", as_dataframe=FALSE, ...)
{
    if (method == "overrepresentation") {

        modules <- assignments(x, as_list=TRUE)

        enrich_res <- reduced_oa(modules, ...)

    } else {
        stop("Enrichment method not recognised")
    }

    if (as_dataframe)  {
        enrich_res <- lapply(enrich_res, function(x) {x@result})
        enrich_res <- do.call("rbind", enrich_res)
    }

    return(enrich_res)
})

setMethod("plotDendro", c("ModularExperiment"),
          function(x, groupLabels = "Module colors", dendroLabels = FALSE,
                       hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                       color_func = WGCNA::labels2colors) {

    colors <- gsub("module_", "", names(assignments(x)))

    is_color <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)),
                     error = function(e) FALSE)
        })
    }

    if (!is.null(color_func) & !all(is_color(colors))) {
        colors <- color_func(colors)
    }

    WGCNA::plotDendroAndColors(dendrogram(x), colors,
                           groupLabels = groupLabels,
                           dendroLabels = dendroLabels, hang = hang,
                           addGuide = addGuide, guideHang = guideHang)
})

setMethod("calcEigengenes", c("ModularExperiment", "matrix"),
          function(x, newdata, ...) {

    eig <- WGCNA::moduleEigengenes(t(newdata), setNames(names(assignments(x)), assignments(x)), ...)
    colnames(eig$eigengenes) <- gsub("ME", "", colnames(eig$eigengenes))
    return(eig$eigengenes)
})

setMethod("calcEigengenes", c("ModularExperiment", "data.frame"), function(x, newdata, ...) {
    return(calcEigengenes(x, as.matrix(newdata), ...))
})

setMethod("calcEigengenes", c("ModularExperiment", "SummarizedExperiment"),
          function(x, newdata, assay_name="normal", ...) {

    eig <- calcEigengenes(x, assay(newdata, assay_name), ...)
    return(.se_to_me(newdata, reduced=as.matrix(eig), assignments=assignments(x)))
})

setMethod("predict", c("ModularExperiment"), function(object, newdata, ...) {
    return(calcEigengenes(object, newdata, ...))
})
