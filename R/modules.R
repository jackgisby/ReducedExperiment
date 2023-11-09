#' Run WGCNA
#' @export
identify_modules <- function(X, powers, min_r_squared = 0.85, max_mean_connectivity = 100,
                             corType = "pearson", networkType = "signed", module_labels = "numbers",
                             ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    wgcna_res <- run_wgcna(assay(X, "normal"), powers=powers, min_r_squared = min_r_squared,
                           corType = corType, networkType=networkType, module_labels = module_labels,
                           seed=seed, ...)
    reduced_set <- ModularExperiment(X, reduced=wgcna_res$E, assignments=wgcna_res$assignments)

    return(reduced_set)
}

.se_to_me <- function(se, assignments, loadings, reduced, center, scale) {
    return(ModularExperiment(assignments=assignments, loadings=loadings,
                             reduced=reduced,
                             assays=assays(se), rowData=rowData(se),
                             colData=colData(se), metadata=metadata(se)))
}

#' Run WGCNA for a data matrix
#' @export
run_wgcna <- function(X, seed=1, powers = 1:40,
                      min_r_squared = 0.85, max_mean_connectivity = 100,
                      corType = "pearson", networkType="signed",
                      module_labels = "numbers", ...) {
    set.seed(seed)

    if (length(powers) > 1) {

        if (corType == "pearson") {
            corFnc <- WGCNA::cor
        } else if (corType == "bicor") {
            corFnc <- WGCNA::bicor
        } else {
            stop("`corType` must be one of 'pearson', 'bicor'")
        }

        threshold <- WGCNA::pickSoftThreshold(t(X), RsquaredCut = min_r_squared, powerVector = powers, corFnc = corFnc, networkType = networkType)

        if (is.null(max_mean_connectivity)) {
            power <- threshold$powerEstimate
        } else {
            power <- min(threshold$Power[threshold$SFT.R.sq > min_r_squared & threshold$mean.k. < max_mean_connectivity])
        }

    } else if (length(powers) == 1) {
        power <- powers
    } else {
        stop("Powers must either be a single integer or a vector of integers to be tested")
    }

    bwms <- WGCNA::blockwiseModules(t(X), power=power, corType=corType, networkType=networkType, ...)

    wgcna_res <- list(
        assignments = bwms$colors,
        E = bwms$MEs
    )

    # TODO: give option for using colours instead `module_labels`
    names(wgcna_res$assignments) <- paste0("factor_", names(wgcna_res$assignments))
    names(wgcna_res$E) <- paste0("factor_", names(wgcna_res$E))

    return(wgcna_res)
}
