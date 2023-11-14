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

    if (!all(sampleNames(object) == rownames(reduced(object))))
        msg <- c(msg, "Reduced data have invalid row names")

    return(if (is.null(msg)) TRUE else msg)
})

setMethod("reduced", "ReducedExperiment", function(x, scale_X=FALSE, center_X=FALSE) {
    return(scale(x@reduced, scale=scale_X, center=center_X))
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

# TODO: setters
setMethod("featureNames", "ReducedExperiment", function(x) {return(names(x))})
setMethod("sampleNames", "ReducedExperiment", function(x) {return(rownames(colData(x)))})

setMethod("show", "ReducedExperiment" ,
          function(object) {
              callNextMethod()
              cat(nComponents(object), " components\n")
          })

setMethod("[", c("ReducedExperiment", "ANY", "ANY", "ANY"),
          function(x, i, j, k, ..., drop=FALSE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    red <- reduced(x)

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

    rowData(x) <- merge(rowData(x), biomart_out, by=gene_id_type, all.x=TRUE)
    stopifnot(all.equal(rownames(x), rownames(rowData(x))))

    return(x)
})
