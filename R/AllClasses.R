
setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("data.frame_OR_NULL", c("data.frame", "NULL"))

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
    representation=representation(loadings="matrix", stability="numeric_OR_NULL", scale="ANY", center="ANY"),
    prototype=prototype(loadings=matrix(), stability=numeric(), scale=TRUE, center=TRUE)
)

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ModularExperiment <- setClass(
    "ModularExperiment",
    contains="ReducedExperiment",
    representation=representation(assignments="character", dendrogram="ANY", threshold="data.frame_OR_NULL"),
    prototype=prototype(assignments=character())
)
