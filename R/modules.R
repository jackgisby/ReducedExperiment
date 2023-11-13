#' Run WGCNA
#' @export
identify_modules <- function(X, ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    wgcna_res <- run_wgcna(assay(X, "normal"), ...)
    reduced_set <- .se_to_me(X, reduced=as.matrix(wgcna_res$E), assignments=wgcna_res$assignments)

    return(reduced_set)
}

.se_to_me <- function(se, assignments, reduced) {
    return(ModularExperiment(assignments=assignments, reduced=reduced,
                             assays=assays(se), rowData=rowData(se),
                             colData=colData(se), metadata=metadata(se)))
}

#' Run WGCNA for a data matrix
#' @export
run_wgcna <- function(X, powers=1:30,
                      min_r_squared=0.85, max_mean_connectivity=100,
                      corType="pearson", networkType="signed",
                      module_labels="numbers", seed=1, verbose = 0, ...) {
    set.seed(seed)

    if (corType == "pearson") {
        cor <- corFnc <- WGCNA::cor
    } else if (corType == "bicor") {
        cor <- corFnc <- WGCNA::bicor
    } else {
        stop("`corType` must be one of 'pearson', 'bicor'")
    }

    if (length(powers) > 1) {

        threshold <- WGCNA::pickSoftThreshold(t(X), RsquaredCut=min_r_squared, powerVector=powers,
                                              corFnc=corFnc, networkType=networkType, verbose=verbose)

        if (is.null(max_mean_connectivity)) {
            power <- threshold$fitIndices$powerEstimate
        } else {
            which_power <- threshold$fitIndices$SFT.R.sq > min_r_squared & threshold$fitIndices$mean.k. < max_mean_connectivity

            if (length(which_power) == 0)
                stop("No power with r_squared > ", min_r_squared, " and mean connectivity < ", max_mean_connectivity)

            power <- min(threshold$fitIndices$Power[which_power])
        }

    } else if (length(powers) == 1) {
        power <- powers
    } else {
        stop("Powers must either be a single integer or a vector of integers to be tested")
    }

    bwms <- WGCNA::blockwiseModules(t(X), power=power, corType=corType, networkType=networkType, verbose=verbose, ...)

    wgcna_res <- list(
        assignments = bwms$colors,
        E = bwms$MEs
    )

    # TODO: give option for using colours instead `module_labels`
    # names(wgcna_res$assignments) <- paste0("factor_", names(wgcna_res$assignments))
    # names(wgcna_res$E) <- paste0("factor_", names(wgcna_res$E))

    return(wgcna_res)
}
