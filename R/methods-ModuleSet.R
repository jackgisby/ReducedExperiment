#' ModularExperiment
#'
#' A class to represent the results of factor analysis
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
ModularExperiment <- function(assignments=character(), ...)
{
    re <- ReducedExperiment(...)
    return(.ModularExperiment(re, assignments=assignments))
}
