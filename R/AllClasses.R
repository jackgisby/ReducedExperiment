#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ReducedExperiment <- setClass(
    "ReducedExperiment",
    contains="SummarizedExperiment",
    representation=representation(reduced="matrix"),
    prototype=prototype(reduced=matrix())
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.FactorisedExperiment <- setClass(
    "FactorisedExperiment",
    contains="ReducedExperiment",
    representation=representation(loadings="matrix", varexp="numeric", center="ANY"),
    prototype=prototype(loadings=matrix(), varexp=numeric(), center=FALSE)
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ModularExperiment <- setClass(
    "ModularExperiment",
    contains="ReducedExperiment",
    representation=representation(assignments="character"),
    prototype=prototype(assignments=character())
)
