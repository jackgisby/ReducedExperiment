#' Run WGCNA
#' @export
identify_modules <- function(X, center_X=TRUE, scale_X=TRUE, assay_name="normal", ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    if (assay_name != "normal") {
        assay(X, "normal") <- assay(X, assay_name)
    }

    if ("transformed" %in% assayNames(X)) warning("Overwriting 'transformed' assay slot in X")
    assay(X, "transformed") <- t(scale(t(assay(X, "normal")), center=center_X, scale=scale_X))

    if (center_X) center_X <- attr(assay(X, "transformed"), "scaled:center")
    if (scale_X) scale_X <- attr(assay(X, "transformed"), "scaled:scale")

    wgcna_res <- run_wgcna(assay(X, "transformed"), return_full_output=FALSE, ...)
    reduced_set <- .se_to_me(X, reduced=wgcna_res$E, loadings=wgcna_res$L, assignments=wgcna_res$assignments, center_X=center_X, scale_X=scale_X, dendrogram=wgcna_res$dendrogram, threshold=wgcna_res$threshold)

    return(reduced_set)
}

.se_to_me <- function(se, reduced, loadings, assignments, center_X, scale_X, dendrogram=NULL, threshold=NULL) {
    return(ModularExperiment(reduced=reduced, loadings=loadings, assignments=assignments,
                             center=center_X, scale=scale_X,
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
                      seed=1, verbose = 0, return_full_output = FALSE,
                      scale_reduced=TRUE, ...) {
    set.seed(seed)

    if (maxBlockSize < nrow(X)) warning("maxBlockSize < total features, module detection will be performed in a block-wise manner")

    if (corType == "pearson") {
        cor <- corFnc <- WGCNA::cor
    } else if (corType == "bicor") {
        cor <- corFnc <- WGCNA::bicor
    } else {
        stop("`corType` must be one of 'pearson', 'bicor'")
    }

    threshold <- suppressWarnings(WGCNA::pickSoftThreshold(
        t(X), RsquaredCut=min_r_squared, powerVector=powers, corFnc=corFnc,
        networkType=networkType, blockSize = maxBlockSize, verbose=verbose)
    )

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

    original_E <- wgcna_res$E

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

    original_E <- wgcna_res$E
    recalculated_E <- .calculate_eigengenes(X, colnames(wgcna_res$E), wgcna_res$assignments, realign = TRUE)
    wgcna_res$E <- recalculated_E$reduced
    wgcna_res$L <- recalculated_E$loadings

    if (scale_reduced) wgcna_res$E <- scale(wgcna_res$E)

    if (return_full_output) {
        return(list("run_wgcna_output" = wgcna_res, "original_E" = original_E, "blockwise_modules_output" = bwms, "pick_soft_threshold_output" = threshold))
    } else {
        return(wgcna_res)
    }
}

.calculate_eigengenes <- function(newdata, module_names, module_assignments, realign = TRUE) {

    red <- data.frame(row.names = colnames(newdata))
    lod <- c()

    for (m in module_names) {

        module_data <- newdata[rownames(newdata) %in% module_assignments[names(module_assignments) == m], ]

        prcomp_res <- prcomp(t(module_data), center = FALSE, scale. = FALSE, rank = 1)

        # Principal components may be anticorrelated, in which case change sign
        align_sign <- ifelse(realign, sign(mean(cor(prcomp_res$x, t(module_data)))), 1)

        red[[m]] <- prcomp_res$x * align_sign

        stopifnot(all(t(module_data) %*% prcomp_res$rotation * align_sign == prcomp_res$x * align_sign))

        lod <- c(lod, setNames(prcomp_res$rotation * align_sign, rownames(module_data)))
    }

    lod <- lod[match(rownames(newdata), names(lod))]

    return(list("reduced" = as.matrix(red), "loadings" = lod))
}

.project_eigengenes <- function(newdata, module_names, module_assignments, lod, min_module_genes) {

    red <- data.frame(row.names = colnames(newdata))

    for (m in module_names) {
        
        module_genes <- module_assignments[names(module_assignments) == m]
        
        if (length(module_genes) < min_module_genes) next
        if (!any(rownames(newdata) %in% module_genes)) next
        
        
        module_data <- newdata[rownames(newdata) %in% module_genes, ]
        module_lod <- lod[names(lod) %in% module_genes]
        stopifnot(all(rownames(module_data) == names(module_lod)))
        
        red[[m]] <- t(module_data) %*% module_lod
    }

    return(as.matrix(red))
}
