#' ReducedExperiment: A container for dimensionally-reduced data
#'
#' @description
#' Inherits from \link[SummarizedExperiment]{SummarizedExperiment}, a
#' container for one or more matrices with features as rows (e.g., genes) and
#' columns as samples. Additional information on features and samples are
#' contained in \link[SummarizedExperiment]{DataFrame} tables. The
#' `ReducedExperiment` extends \link[SummarizedExperiment]{SummarizedExperiment}
#' by additionally providing access to a "reduced" data matrix, in which rows
#' represent samples and columns represent a second set of dimensionally-reduced
#' features.
#'
#' The methods available for \link[SummarizedExperiment]{SummarizedExperiment}
#' objects are also available for `ReducedExperiment` and its children, including
#' \link{FactorisedExperiment} and \link{ModularExperiment}
#'
#' Typically, `ReducedExperiment` objects contain two main assays. The first is,
#' by default, named "normal" and contains some type of normalised data,
#' such as gene expression data. The second is "transformed", which is typically
#' the result of applying scaling and/or centering to the normalised data
#' matrix.
#'
#' @param reduced A data matrix, usually the result of some type of
#' dimensionality-reduction, with rows representing samples and columns
#' representing a new set of features.
#'
#' @param scale Either a boolean, representing whether or not the original data
#' has been scaled to unit variance, or a numeric vector indicating the
#' standard deviations of the original features (as produced by
#' \link[base]{scale}.)
#'
#' @param center Either a boolean, representing whether or not the original data
#' has been centered to have a mean of 0, or a numeric vector indicating the
#' means of the original features (as produced by
#' \link[base]{scale}.)
#'
#' @import SummarizedExperiment
#'
#' @export
ReducedExperiment <- function(
        reduced = new("matrix"),
        scale = TRUE,
        center = TRUE,
        ...)
{
    se <- SummarizedExperiment::SummarizedExperiment(...)
    return(.ReducedExperiment(se, reduced=reduced, scale=scale, center=center))
}

S4Vectors::setValidity2("ReducedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Check sampleNames method matches with other methods
    if (!identical(sampleNames(object), colnames(object)))
        msg <- c(msg, "Reduced data have invalid row names")

    # Check featureNames matches with other methods
    if (!identical(featureNames(object), names(object)))
        msg <- c(msg, "Feature names do not match with names")
    if (!identical(featureNames(object), rownames(object)))
        msg <- c(msg, "Feature names do not match with rownames")

    # Center / Scale - check names matches feature names
    if (!is.logical(object@scale)) {
        if (!identical(rownames(object), names(object@scale)))
            msg <- c(msg, "Scaling vector has invalid names")
    }
    if (!is.logical(object@center)) {
        if (!identical(rownames(object), names(object@center)))
            msg <- c(msg, "Centering vector has invalid names")
    }

    # Check reduced matrix
    if (obj_dims[2] != dim(reduced(object))[1])
        msg <- c(msg, "Reduced data have invalid row dimensions")

    if (!identical(sampleNames(object), rownames(reduced(object))))
        msg <- c(msg, "Reduced data have invalid row names")

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("reduced", "ReducedExperiment", function(x, scale_reduced=FALSE, center_reduced=FALSE) {
    return(scale(x@reduced, scale=scale_reduced, center=center_reduced))
})

setReplaceMethod("reduced", "ReducedExperiment", function(x, value) {
    x@reduced <- value
    validObject(x)
    return(x)
})

setMethod("componentNames", "ReducedExperiment", function(x) {return(colnames(x@reduced))})

setReplaceMethod("componentNames", "ReducedExperiment", function(x, value) {
    colnames(x@reduced) <- value
    validObject(x)
    return(x)
})

setMethod("featureNames", "ReducedExperiment", function(x) {return(names(x))})

setReplaceMethod("names", "ReducedExperiment", function(x, value) {
    x <- callNextMethod(x, value)
    if (!is.logical(x@scale)) names(x@scale) <- value
    if (!is.logical(x@center)) names(x@center) <- value
    validObject(x)
    return(x)
})
setReplaceMethod("rownames", "ReducedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})
setReplaceMethod("featureNames", "ReducedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})


setMethod("sampleNames", "ReducedExperiment", function(x) {return(colnames(x))})

setReplaceMethod("sampleNames", "ReducedExperiment", function(x, value) {
    rownames(x@reduced) <- colnames(x) <- value
    validObject(x)
    return(x)
})

setMethod("show", "ReducedExperiment" ,
          function(object) {
              callNextMethod()
              cat(nComponents(object), "components\n")
          })

setMethod("[", c("ReducedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    red <- x@reduced
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
        if (!is.logical(center)) center <- center[i,drop=FALSE]
        if (!is.logical(scale)) scale <- scale[i,drop=FALSE]
    }

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
                k, componentNames(x), fmt
            )
        }

        k <- as.vector(k)
        red <- red[,k,drop=FALSE]
    }

    out <- callNextMethod(x, i, j, ...)
    BiocGenerics:::replaceSlots(out, reduced=red, center=center, scale=scale, check=FALSE)
})

setMethod("dim", "ReducedExperiment", function(x) {
    out <- c(callNextMethod(x), ncol(reduced(x)))
    names(out) <- c("Features", "Samples", "Components")
    return(out)
})

setMethod("nComponents", "ReducedExperiment", function(x) {dim(x)[3]})
setMethod("nSamples", "ReducedExperiment", function(x) {dim(x)[2]})
setMethod("nFeatures", "ReducedExperiment", function(x) {dim(x)[1]})

#' @import biomaRt
setMethod("getGeneIDs", "ReducedExperiment", function(
    x,
    gene_id_col="rownames",
    gene_id_type="ensembl_gene_id",
    ids_to_get=c("hgnc_symbol", "entrezgene_id"),
    dataset="hsapiens_gene_ensembl",
    mart=NULL)
{
    if (gene_id_col == "rownames") {
        rowData(x)[[gene_id_type]] <- rownames(x)
    } else {
        rowData(x)[[gene_id_type]] <- rowData(x)[[gene_id_col]]
    }

    gene_ids <- rowData(x)[[gene_id_type]]

    if (is.null(mart)) {
        mart <- biomaRt::useEnsembl(biomart="genes", dataset=dataset)
    }

    biomart_out <- biomaRt::getBM(filters = gene_id_type,
                                  attributes = c(gene_id_type, ids_to_get),
                                  values = gene_ids, mart = mart)

    biomart_out <- biomart_out[which(!duplicated(biomart_out[[gene_id_type]])),]
    rownames(biomart_out) <- biomart_out[[gene_id_type]]

    # TODO: This approach can probably be improved
    row_data_merged <- merge(rowData(x), biomart_out, by=gene_id_type, all.x=TRUE)
    rownames(row_data_merged) <- row_data_merged[[gene_id_type]]
    row_data_merged <- row_data_merged[match(rowData(x)[[gene_id_type]], row_data_merged[[gene_id_type]]) ,]

    stopifnot(identical(row_data_merged[[gene_id_type]], rownames(row_data_merged)))
    stopifnot(identical(rownames(x), rownames(row_data_merged)))

    rowData(x) <- row_data_merged

    return(x)
})
