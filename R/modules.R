#' Run WGCNA
#' @export
run_wgcna <- function(X, ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list(normal = X))
    }

    wgcna_res <- run_wgcna(X, ...)
    reduced_set <- ModularExperiment(X, reduced=wgcna_res$E,
                                     assignments=wgcna_res$assignments)

    return(reduced_set)
}
