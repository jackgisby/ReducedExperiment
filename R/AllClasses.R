setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("data.frame_OR_NULL", c("data.frame", "NULL"))
setClassUnion("logical_OR_numeric", c("logical", "numeric"))

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @rdname reduced_experiment
#' @exportClass ReducedExperiment
.ReducedExperiment <- setClass(
    "ReducedExperiment",
    contains = "SummarizedExperiment",
    representation = representation(
        reduced = "matrix",
        scale = "logical_OR_numeric",
        center = "logical_OR_numeric"
    ),
    prototype = prototype(
        reduced = matrix(),
        scale = TRUE,
        center = TRUE
    )
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @rdname factorised_experiment
#' @exportClass FactorisedExperiment
.FactorisedExperiment <- setClass(
    "FactorisedExperiment",
    contains = "ReducedExperiment",
    representation = representation(
        loadings = "matrix",
        stability = "numeric_OR_NULL"
    ),
    prototype = prototype(
        loadings = matrix(),
        stability = numeric()
    )
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @rdname modular_experiment
#' @exportClass ModularExperiment
.ModularExperiment <- setClass(
    "ModularExperiment",
    contains = "ReducedExperiment",
    representation = representation(
        loadings = "numeric_OR_NULL",
        assignments = "character",
        dendrogram = "ANY",
        threshold = "data.frame_OR_NULL"
    ),
    prototype = prototype(
        loadings = numeric(),
        assignments = character()
    )
)
