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
    representation=representation(loadings="matrix", varexp="numeric", scale="ANY", center="ANY"),
    prototype=prototype(loadings=matrix(), varexp=numeric(), scale=TRUE, center=TRUE)
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ModularExperiment <- setClass(
    "ModularExperiment",
    contains="ReducedExperiment",
    representation=representation(assignments="character"),
    prototype=prototype(assignments=character())
)
