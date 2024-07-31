#' ReducedExperiment: A container for dimensionally-reduced data
#'
#' @description
#' Inherits from \link[SummarizedExperiment]{SummarizedExperiment}, a
#' container for one or more matrices with features as rows (e.g., genes) and
#' columns as samples. Additional information on features and samples are
#' contained in \link[S4Vectors]{DataFrame} tables. The
#' `ReducedExperiment` extends \link[SummarizedExperiment]{SummarizedExperiment}
#' by additionally providing access to a "reduced" data matrix, in which rows
#' represent samples and columns represent a second set of dimensionally-reduced
#' features.
#'
#' The methods available for \link[SummarizedExperiment]{SummarizedExperiment}
#' objects are also available for `ReducedExperiment` and its children, which
#' include \link[ReducedExperiment]{FactorisedExperiment} and
#' \link[ReducedExperiment]{ModularExperiment}.
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
#' @param ... Additional arguments to be passed to
#' \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @returns Constructor method returns a
#' \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @seealso [ReducedExperiment::FactorisedExperiment()],
#' [ReducedExperiment::ModularExperiment()]
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' # In this case we use random assay and reduced data, but in
#' # practice these will likely be the result of applying some kind of
#' # dimensionality-reduction method to the assay data (e.g., gene
#' # expression data) from some study.
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' re
#'
#' @import SummarizedExperiment
#'
#' @rdname reduced_experiment
#' @export
ReducedExperiment <- function(
        reduced = new("matrix"),
        scale = TRUE,
        center = TRUE,
        ...) {
    se <- SummarizedExperiment::SummarizedExperiment(...)

    return(.ReducedExperiment(
        se,
        reduced = reduced,
        scale = scale,
        center = center
    ))
}

S4Vectors::setValidity2("ReducedExperiment", function(object) {
    msg <- NULL

    obj_dims <- dim(object)

    # Check sampleNames method matches with other methods
    if (!identical(sampleNames(object), colnames(object))) {
        msg <- c(msg, "Reduced data have invalid row names")
    }

    # Check featureNames matches with other methods
    if (!identical(featureNames(object), names(object))) {
        msg <- c(msg, "Feature names do not match with names")
    }
    if (!identical(featureNames(object), rownames(object))) {
        msg <- c(msg, "Feature names do not match with rownames")
    }

    # Center / Scale - check names matches feature names
    if (!is.logical(object@scale)) {
        if (!identical(rownames(object), names(object@scale))) {
            msg <- c(msg, "Scaling vector has invalid names")
        }
    }
    if (!is.logical(object@center)) {
        if (!identical(rownames(object), names(object@center))) {
            msg <- c(msg, "Centering vector has invalid names")
        }
    }

    # Check reduced matrix
    if (obj_dims[2] != dim(reduced(object))[1]) {
        msg <- c(msg, "Reduced data have invalid row dimensions")
    }

    if (!identical(sampleNames(object), rownames(reduced(object)))) {
        msg <- c(msg, "Reduced data have invalid row names")
    }

    return(if (is.null(msg)) TRUE else msg)
})

#' Get and set reduced data
#'
#' Retrieves the reduced data matrix, with samples as rows and reduced
#' components as columns. Setter method can be used to replace or modify
#' all or part of the matrix.
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object or its
#' children.
#'
#' @param scale_reduced If `TRUE`, data will be scaled column-wise to have a
#' standard deviation of 0.
#'
#' @param center_reduced If `TRUE`, data will be centered column-wise to have a
#' mean of 0.
#'
#' @param value New value to replace existing reduced data matrix.
#'
#' @returns A matrix with samples/observations as rows and columns referring
#' to the dimensionally-reduced components.
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' stopifnot(all.equal(reduced(re), rand_reduced_data))
#'
#' print(paste0("Reduced data at position (2, 2): ", reduced(re)[2, 2]))
#' reduced(re)[2, 2] <- 5
#' print(paste0("Reduced data at position (2, 2): ", reduced(re)[2, 2]))
#'
#' @seealso [ReducedExperiment::ReducedExperiment()]
#'
#' @rdname reduced
#' @name reduced
#' @aliases reduced<-
#' @export reduced
NULL

#' @rdname reduced
#' @export
setMethod("reduced", "ReducedExperiment", function(
        object,
        scale_reduced = FALSE,
        center_reduced = FALSE) {
    return(scale(
        object@reduced,
        scale = scale_reduced,
        center = center_reduced
    ))
})

#' @rdname reduced
#' @export
setReplaceMethod("reduced", "ReducedExperiment", function(object, value) {
    object@reduced <- value
    validObject(object)
    return(object)
})

#' Get names of dimensionally-reduced components
#'
#' Retrieves the feature names post-dimensionality reduction In the case of
#' module analysis, these are the names of the gene modules; in the case of
#' factor analysis, these are the names of the factors.
#'
#' @param object A \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param value New value to replace existing names.
#'
#' @returns A vector containing the names of the components.
#'
#' @details
#' `componentNames` is valid for all \link[ReducedExperiment]{ReducedExperiment}
#' objects, whereas `moduleNames` is only valid for
#' \link[ReducedExperiment]{ModularExperiment}s.
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' stopifnot(all.equal(componentNames(re), colnames(rand_reduced_data)))
#'
#' print(paste0("Component name at [2]: ", componentNames(re)[2]))
#' componentNames(re)[2] <- "custom_component_name"
#' print(paste0("Component name at [2]: ", componentNames(re)[2]))
#'
#' @rdname component_names
#' @name componentNames
#' @aliases componentNames<- moduleNames moduleNames<-
#' @export componentNames
NULL

#' @rdname component_names
#' @export
setMethod("componentNames", "ReducedExperiment", function(object) {
    return(colnames(object@reduced))
})

#' @rdname component_names
#' @export
setReplaceMethod("componentNames", "ReducedExperiment", function(
        object,
        value) {
    colnames(object@reduced) <- value
    validObject(object)
    return(object)
})

#' Get feature names
#'
#' Retrieves feature names (rownames, usually genes).
#'
#' @param x \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param value New value to replace existing names.
#'
#' @returns A vector containing the names of the features.
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' stopifnot(all.equal(featureNames(re), rownames(rand_assay_data)))
#' stopifnot(all.equal(rownames(re), rownames(rand_assay_data)))
#' stopifnot(all.equal(names(re), rownames(rand_assay_data)))
#'
#' print(paste0("Feature name at position 55: ", featureNames(re)[55]))
#' featureNames(re)[55] <- "custom_feature_name"
#' print(paste0("Reduced data at position 55: ", featureNames(re)[55]))
#'
#' @rdname feature_names
#' @name featureNames
#' @aliases featureNames<-
#' @export featureNames
NULL

#' @rdname feature_names
#' @export
setMethod("featureNames", "ReducedExperiment", function(x) {
    return(names(x))
})

#' @rdname feature_names
#' @export
setReplaceMethod("names", "ReducedExperiment", function(x, value) {
    x <- callNextMethod(x, value)
    if (!is.logical(x@scale)) names(x@scale) <- value
    if (!is.logical(x@center)) names(x@center) <- value
    validObject(x)
    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("rownames", "ReducedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' @rdname feature_names
#' @export
setReplaceMethod("featureNames", "ReducedExperiment", function(x, value) {
    names(x) <- value
    return(x)
})

#' Get sample names
#'
#' Retrieves sample names (colnames).
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param value New value to replace existing names.
#' @returns A vector containing the names of the features.
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' stopifnot(all.equal(sampleNames(re), colnames(rand_assay_data)))
#' stopifnot(all.equal(colnames(re), colnames(rand_assay_data)))
#'
#' print(paste0("Sample name at [80]: ", sampleNames(re)[80]))
#' sampleNames(re)[80] <- "custom_feature_name"
#' print(paste0("Sample data at [80]: ", sampleNames(re)[80]))
#'
#' @rdname sample_names
#' @name sampleNames
#' @aliases sampleNames<-
#' @export sampleNames
NULL

#' @rdname sample_names
#' @export
setMethod("sampleNames", "ReducedExperiment", function(object) {
    return(colnames(object))
})

#' @rdname sample_names
#' @export
setReplaceMethod("sampleNames", "ReducedExperiment", function(object, value) {
    rownames(object@reduced) <- colnames(object) <- value
    validObject(object)
    return(object)
})

#' Prints a summary of a ReducedExperiment object
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @returns A character summary describing the object.
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of factors
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' # Equivalent to `show(re)`
#' re
#'
#' @rdname show
#' @name show
NULL

#' @rdname show
#' @export
setMethod(
    "show", "ReducedExperiment",
    function(object) {
        callNextMethod()
        cat(nComponents(object), "components\n")
    }
)

#' Command line completion for `$`
#'
#' @description
#' Command line completion for `$`.
#' This function is not intended to be used directly by users but provides
#' auto-completion capabilities.
#' Autocompletes based on column data names (i.e., the column names of the
#' `colData`).
#'
#' @param x The \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param pattern Search pattern.
#'
#' @return The names of the matching columns of `colData`.
#'
#' @seealso [utils::.DollarNames()]
#'
#' @importFrom utils .DollarNames
#'
#' @rdname dollar_names
#' @name dollar_names
NULL

#' @rdname dollar_names
#' @export
.DollarNames.ReducedExperiment <- function(x, pattern = "") {
    grep(pattern, colnames(colData(x)), value = TRUE)
}

#' Extract and replace parts of ReducedExperiment objects
#'
#' @description
#' Method permits slicing of \link[ReducedExperiment]{ReducedExperiment}
#' objects.
#'
#' @param x \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param i Slicing by rows (features, usually genes).
#'
#' @param j Slicing by columns (samples/observations).
#'
#' @param k Slicing by reduced dimensions.
#'
#' @param drop Included for consistency with other slicing methods.
#'
#' @param ... Additional arguments to be passed to the parent method.
#'
#' @returns A \link[ReducedExperiment]{ReducedExperiment} object sliced by
#' rows (`i`), columns (`j`) and components (`k`).
#'
#' @examples
#' # Create randomised data with the following dimensions
#' i <- 300 # Number of features
#' j <- 100 # Number of samples
#' k <- 10 # Number of components (i.e., factors/modules)
#'
#' rand_assay_data <- ReducedExperiment:::.makeRandomData(i, j, "gene", "sample")
#' rand_reduced_data <- ReducedExperiment:::.makeRandomData(j, k, "sample", "component")
#'
#' # Create a randomised ReducedExperiment container
#' re <- ReducedExperiment(
#'     assays = list("normal" = rand_assay_data),
#'     reduced = rand_reduced_data
#' )
#'
#' # Slice our object by rows (1:50), columns (1:20) and components (1:5)
#' # re[i, j, k, ...]
#' re[1:50, 1:20, 1:5]
#'
#' @rdname slice
#' @name slice
NULL

#' @rdname slice
#' @export
setMethod(
    "[", c("ReducedExperiment", "ANY", "ANY", "ANY"),
    function(x, i, j, k, ..., drop = FALSE) {
        object <- x

        if (1L != length(drop) || (!missing(drop) && drop)) {
            warning("'drop' ignored '[,", class(object), ",ANY,ANY-method'")
        }

        red <- object@reduced
        center <- object@center
        scale <- object@scale

        if (!missing(i)) {
            if (is.character(i)) {
                fmt <- paste0(
                    "<", class(object),
                    ">[i,] index out of bounds: %s"
                )
                i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    i, rownames(object), fmt
                )
            }
            i <- as.vector(i)
            if (!is.logical(center)) center <- center[i, drop = FALSE]
            if (!is.logical(scale)) scale <- scale[i, drop = FALSE]
        }

        if (!missing(j)) {
            if (is.character(j)) {
                fmt <- paste0(
                    "<", class(object),
                    ">[,j] index out of bounds: %s"
                )
                j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    j, colnames(object), fmt
                )
            }
            j <- as.vector(j)
            red <- red[j, , drop = FALSE]
        }

        if (!missing(k)) {
            if (is.character(k)) {
                fmt <- paste0(
                    "<", class(object),
                    ">[k,] index out of bounds: %s"
                )
                k <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                    k, componentNames(object), fmt
                )
            }

            k <- as.vector(k)
            red <- red[, k, drop = FALSE]
        }

        out <- callNextMethod(object, i, j, ...)
        BiocGenerics:::replaceSlots(
            out,
            reduced = red,
            center = center,
            scale = scale,
            check = FALSE
        )
    }
)

#' Combine ReducedExperiment objects by columns
#'
#' @description
#' Combines \link[ReducedExperiment]{ReducedExperiment} objects by columns
#' (samples). This method assumes that objects have identical features and
#' components (i.e., factors or modules). If they are not, an error is returned.
#'
#' So, this means that the feature-level slots should be equivalent, for example
#' the assay rownames and the rownames of the factor `loadings` available in
#' \link[ReducedExperiment]{FactorisedExperiment} objects. The component slots
#' should also be equivalent, such as the column names of the `reduced` matrix
#' or the column names of the aformentioned factor `loadings` matrix.
#'
#' @param ... A series of \link[ReducedExperiment]{ReducedExperiment} objects
#' to be combined. See
#' \link[SummarizedExperiment]{cbind,SummarizedExperiment-method} for
#' further details.
#'
#' @param deparse.level Integer, see \link[base]{cbind}.
#'
#' @returns Returns a single \link[ReducedExperiment]{ReducedExperiment} object
#' containing all of the columns in the objects passed to `cbind`.
#'
#' @details
#' The \link[SummarizedExperiment]{SummarizedExperiment} package includes
#' separate methods for `cbind`
#' (\link[SummarizedExperiment]{cbind,SummarizedExperiment-method}) and
#' (\link[SummarizedExperiment]{combineRows}). The latter is supposed to be
#' more flexible, permitting differences in the number and identity of the rows.
#' For \link[ReducedExperiment]{ReducedExperiment} objects we only implement a
#' single, less flexible, method that assumes the rows and components
#' (i.e., factors or modules) are identical across objects. Attempting to apply
#' `combineRows` to a \link[ReducedExperiment]{ReducedExperiment} object will
#' result in the objects being treated as if they were
#' \link[SummarizedExperiment]{SummarizedExperiment}s, and a single
#' \link[SummarizedExperiment]{SummarizedExperiment} object will be returned.
#'
#' @examples
#' # Create randomised containers with different numbers of samples
#' i <- 300 # Number of features
#' k <- 10 # Number of components (i.e., factors/modules)
#'
#' # Same features and components, different samples (30 vs. 50 columns)
#' re_1 <- ReducedExperiment:::.createRandomisedReducedExperiment(i, 50, k)
#' re_2 <- ReducedExperiment:::.createRandomisedReducedExperiment(i, 30, k)
#'
#' # Make a new object with 80 columns
#' cbind(re_1, re_2)
#'
#' # We can apply combineRows to `ReducedExperiment` objects but the resulting
#' # object will be a `SummarizedExperiment`
#' combineRows(re_1, re_2)
#'
#' @seealso [base::cbind()], \link[SummarizedExperiment]{cbind,SummarizedExperiment-method}
#'
#' @rdname cbind
#' @name cbind
NULL

#' @rdname cbind
#' @export
setMethod("cbind", "ReducedExperiment", function(..., deparse.level = 1) {
    args <- list(...)

    reduced <- do.call(rbind, lapply(args, reduced))

    std_slots_equal <- vapply(args, function(re) {
        return(identical(re@scale, args[[1]]@scale) &
            identical(re@center, args[[1]]@center))
    }, FUN.VALUE = FALSE)

    if (all(std_slots_equal)) {
        args[[1]] <- BiocGenerics:::replaceSlots(
            args[[1]],
            reduced = reduced,
            check = FALSE
        )
    } else {
        stop("Row bind expects scale and center slots are equal")
    }

    args[["deparse.level"]] <- deparse.level

    return(do.call(callNextMethod, args))
})

#' Get the dimensions of a Reducedexperiment object
#'
#' @param x \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @returns Returns a named vector containing the dimensions of the samples,
#' features and reduced dimensions.
#'
#' @examples
#' # Create a randomised ReducedExperiment
#' re <- ReducedExperiment:::.createRandomisedReducedExperiment(100, 50, 10)
#'
#' # Get the dimensions
#' dim(re)
#'
#' @export
setMethod("dim", "ReducedExperiment", function(x) {
    object <- x

    out <- c(callNextMethod(object), ncol(object@reduced))
    names(out) <- c("Features", "Samples", "Components")
    return(out)
})

#' Prints individual lengths of samples, components and features
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @returns The number of samples (`nSamples`), features (`nFeatures`)
#' or dimensionally-reduced components (`nComponents`) are returned.
#'
#' @examples
#' # Create a randomised ReducedExperiment
#' re <- ReducedExperiment:::.createRandomisedReducedExperiment(100, 50, 10)
#'
#' # Get the dimensions
#' nComponents(re) # 10
#' nSamples(re) # 50
#' nFeatures(re) # 10
#'
#' # For a ModularExperiment we can alternatively use nModules
#' me <- ReducedExperiment:::.createRandomisedModularExperiment(100, 50, 10)
#' nComponents(me) # 10
#' nModules(me) # 10
#'
#' @seealso \link[ReducedExperiment]{dim,ReducedExperiment-method}
#'
#' @rdname individual_dims
#' @name individual_dim
#' @aliases nComponents nFeatures nModules nSamples
NULL

#' @rdname individual_dims
#' @export
setMethod("nComponents", "ReducedExperiment", function(object) {
    dim(object)[3]
})

#' @rdname individual_dims
#' @export
setMethod("nSamples", "ReducedExperiment", function(object) {
    dim(object)[2]
})

#' @rdname individual_dims
#' @export
setMethod("nFeatures", "ReducedExperiment", function(object) {
    dim(object)[1]
})

#' Gets alternative gene annotations from biomaRt
#'
#' Uses \link[biomaRt]{getBM} to get alternative gene IDs for
#' \link[ReducedExperiment]{ReducedExperiment} objects. The new annotations
#' are added as columns to the input object's `rowData`
#'
#' @param object \link[ReducedExperiment]{ReducedExperiment} object.
#'
#' @param gene_id_col The column in `rowData(object)` that will be used to
#' query biomaRt. Setting this to "rownames" instead uses `rownames(object)`
#' for matching.
#'
#' @param gene_id_type The type of attribute to be used to query with biomaRt.
#' See the `filters` argument of \link[biomaRt]{getBM}.
#'
#' @param ids_to_get The type of attribute to get from biomaRt.
#' See the `attributes` argument of \link[biomaRt]{getBM}.
#'
#' @param dataset The Ensembl dataset to retrieve. See the `dataset` argument
#' of \link[biomaRt]{useEnsembl}. If `mart` is not NULL, this argument is
#' ignored.
#'
#' @param mart An optional mart object to use. See the `mart` argument of
#' \link[biomaRt]{getBM}. If provided, this object is used to query biomart
#' for the conversion of gene IDs. If `biomart_out` is not NULL, this argument
#' is ignored.
#'
#' @param biomart_out An optional `data.frame` containing the output of a call
#' to \link[biomaRt]{getBM}. If provided, this object is used for the conversion
#' of gene IDs.
#'
#' @returns Returns the original object, with additional variables added to
#' the `rowData` slot.
#'
#' @seealso [biomaRt::useEnsembl()], [biomaRt::getBM()]
#'
#' @examples
#' set.seed(2)
#' airway <- ReducedExperiment:::.get_airway_data(n_features = 500)
#'
#' set.seed(1)
#' airway_fe <- estimate_factors(airway, nc = 2, use_stability = FALSE, method = "imax")
#'
#' # rowData before getting additional gene IDs
#' rowData(airway_fe)
#'
#' # For this example we run `getGeneIDs` using a preloaded biomart query
#' # (`biomart_out`) to avoid actually querying ensembl during testing
#' # Note: do not use this file for your actual data
#' biomart_out <- read.csv(system.file(
#'     "extdata",
#'     "biomart_out.csv",
#'     package = "ReducedExperiment"
#' ))
#' airway_fe <- getGeneIDs(airway_fe, biomart_out = biomart_out)
#'
#' # rowData after getting additional gene IDs
#' rowData(airway_fe)
#'
#' @import biomaRt
#'
#' @rdname get_gene_ids
#' @name getGeneIDs
#' @export getGeneIDs
NULL

#' @rdname get_gene_ids
#' @export
setMethod("getGeneIDs", "ReducedExperiment", function(
        object,
        gene_id_col = "rownames",
        gene_id_type = "ensembl_gene_id",
        ids_to_get = c("hgnc_symbol", "entrezgene_id"),
        dataset = "hsapiens_gene_ensembl",
        mart = NULL,
        biomart_out = NULL) {
    if (gene_id_col == "rownames") {
        rowData(object)[[gene_id_type]] <- rownames(object)
    } else {
        rowData(object)[[gene_id_type]] <- rowData(object)[[gene_id_col]]
    }

    gene_ids <- rowData(object)[[gene_id_type]]

    if (is.null(biomart_out)) {
        if (is.null(mart)) {
            mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
        }

        biomart_out <- biomaRt::getBM(
            filters = gene_id_type,
            attributes = c(gene_id_type, ids_to_get),
            values = gene_ids, mart = mart
        )
    }

    biomart_out <-
        biomart_out[which(!duplicated(biomart_out[[gene_id_type]])), ]
    rownames(biomart_out) <- biomart_out[[gene_id_type]]

    # TODO: This approach can probably be improved
    row_data_merged <- merge(rowData(object), biomart_out,
        by = gene_id_type, all.x = TRUE
    )
    rownames(row_data_merged) <- row_data_merged[[gene_id_type]]
    row_data_merged <- row_data_merged[match(
        rowData(object)[[gene_id_type]],
        row_data_merged[[gene_id_type]]
    ), ]

    stopifnot(identical(
        row_data_merged[[gene_id_type]],
        rownames(row_data_merged)
    ))
    stopifnot(identical(rownames(object), rownames(row_data_merged)))

    rowData(object) <- row_data_merged

    return(object)
})
