#' Run WGCNA
#' @export
identify_modules <- function(X, assay_name="normal", ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    if (assay_name != "normal") {
        assay(X, "normal") <- assay(X, assay_name)
    }

    wgcna_res <- run_wgcna(assay(X, "normal"), return_full_output=FALSE, ...)
    reduced_set <- .se_to_me(X, reduced=as.matrix(wgcna_res$E), assignments=wgcna_res$assignments, dendrogram=wgcna_res$dendrogram, threshold=wgcna_res$threshold)

    return(reduced_set)
}

.se_to_me <- function(se, assignments, reduced, dendrogram=NULL, threshold=NULL) {
    return(ModularExperiment(assignments=assignments, reduced=reduced,
                             dendrogram=dendrogram, threshold=threshold,
                             assays=assays(se), rowData=rowData(se),
                             colData=colData(se), metadata=metadata(se)))
}

#' Run WGCNA for a data matrix
#' @export
run_wgcna <- function(X, powers=1:30,
                      min_r_squared=0.85, max_mean_connectivity=100,
                      corType="pearson", networkType="signed",
                      module_labels="numbers", maxBlockSize = 30000,
                      seed=1, verbose = 0, return_full_output = FALSE, ...) {
    set.seed(seed)

    if (maxBlockSize < nrow(X)) warning("maxBlockSize < total features, module detection will be performed in a block-wise manner")

    if (corType == "pearson") {
        cor <- corFnc <- WGCNA::cor
    } else if (corType == "bicor") {
        cor <- corFnc <- WGCNA::bicor
    } else {
        stop("`corType` must be one of 'pearson', 'bicor'")
    }

    threshold <- WGCNA::pickSoftThreshold(t(X), RsquaredCut=min_r_squared, powerVector=powers,
                                          corFnc=corFnc, networkType=networkType,
                                          blockSize = maxBlockSize, verbose=verbose)

    if (length(powers) > 1) {

        if (is.null(max_mean_connectivity)) {
            power <- threshold$fitIndices$powerEstimate
        } else {
            which_power <- which(threshold$fitIndices$SFT.R.sq > min_r_squared & threshold$fitIndices$mean.k. < max_mean_connectivity)

            if (length(which_power) == 0)
                stop("No power with r_squared > ", min_r_squared, " and mean connectivity < ", max_mean_connectivity)

            power <- min(threshold$fitIndices$Power[which_power])
        }

    } else if (length(powers) == 1) {
        power <- powers
    }

    threshold$fitIndices$selected_power <- power

    bwms <- WGCNA::blockwiseModules(t(X), power=power, corType=corType, networkType=networkType,
                                    maxBlockSize=maxBlockSize, verbose=verbose, ...)

    wgcna_res <- list(
        assignments = bwms$colors,
        E = bwms$MEs,
        dendrogram = if (length(bwms$dendrograms) == 1) bwms$dendrograms[[1]] else NULL,
        threshold = threshold$fitIndices
    )

    colnames(wgcna_res$E) <- gsub("ME", "", colnames(wgcna_res$E))

    if (module_labels == "numbers") {

        .colors2numbers <- function(colors) {

            color_table <- table(colors)
            color_table <- color_table[order(color_table, decreasing = TRUE)]
            color_table <- color_table[which(names(color_table) != "grey")]

            color_table <- setNames(1:length(color_table), names(color_table))
            color_table <- c(color_table, "grey" = 0)

            return(color_table)
        }

        converter <- .colors2numbers(wgcna_res$assignments)
        for (i in 1:length(wgcna_res$assignments)) {
            wgcna_res$assignments[i] <- converter[which(names(converter) == wgcna_res$assignments[i])]
        }
        for (i in 1:length(colnames(wgcna_res$E))) {
            colnames(wgcna_res$E)[i] <- converter[which(names(converter) == colnames(wgcna_res$E)[i])]
        }

    } else if (module_labels != "colors" & module_labels != "colours") {
        stop("Value of `module_labels` does not correspond to a valid option")
    }

    wgcna_res$assignments <- setNames(names(wgcna_res$assignments), paste0("module_", wgcna_res$assignments))
    colnames(wgcna_res$E) <- paste0("module_", colnames(wgcna_res$E))
    wgcna_res$E <- wgcna_res$E[, order(colnames(wgcna_res$E))]

    if (return_full_output) {
        return(list("run_wgcna_output" = wgcna_res, "blockwise_modules_output" = bwms, "pick_soft_threshold_output" = threshold))
    } else {
        return(wgcna_res)
    }
}
