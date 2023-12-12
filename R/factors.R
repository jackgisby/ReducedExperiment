#' Apply dimensionality reduction using ICA
#' @export
estimate_factors <- function(X, nc, center_X=TRUE, scale_X=FALSE, ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    if ("transformed" %in% assayNames(X)) warning("Overwriting 'transformed' assay slot in X")
    assay(X, "transformed") <- t(scale(t(assay(X, "normal")), center=center_X, scale=scale_X))

    if (center_X) center_X <- attr(assay(X, "transformed"), "scaled:center")
    if (scale_X) scale_X <- attr(assay(X, "transformed"), "scaled:scale")

    ica_res <- run_ica(assay(X, "transformed"), nc=nc,
                       center_X=FALSE, scale_X=FALSE, ...)

    return(.se_to_fe(X, reduced=ica_res$M, loadings=ica_res$S, stability=ica_res$stab, center=center_X, scale=scale_X))
}

.se_to_fe <- function(se, reduced, loadings, stability, center, scale) {
    return(FactorisedExperiment(loadings=loadings, stability=stability,
                                center=center, scale=scale, reduced=reduced,
                                assays=assays(se), rowData=rowData(se),
                                colData=colData(se), metadata=metadata(se)))
}

#' Run ICA for a data matrix
#' @export
run_ica <- function(X, nc, use_stability=TRUE, resample=FALSE,
                    method=ifelse(use_stability, "fast", "imax"),
                    center_X=TRUE, scale_X=FALSE,
                    reorient_skewed=TRUE, seed=1,
                    scale_components=TRUE, BPPARAM = SerialParam(),
                    ...) {
    set.seed(seed)

    if (center_X | scale_X)
        {X <- t(scale(t(X), center=center_X, scale=scale_X))}

    if (use_stability) {
        ica_res <- .stability_ica(X, nc=nc, resample=resample, method=method, BPPARAM=BPPARAM, ...)
    } else {
        if (resample) stop("Cannot use resampling approach when `use_stability` is FALSE")
        ica_res <- list(S = ica::ica(X, nc=nc, method=method, center=FALSE, ...)$S)
    }

    # Reorient and scale factors before recalculating M
    if (reorient_skewed) {ica_res$S <- .reorient_factors(ica_res$S)}
    if (scale_components) ica_res$S <- scale(ica_res$S, center = FALSE)
    ica_res$M <- .project_ica(X, ica_res$S)

    # Add factors / sample names
    rownames(ica_res$M) <- colnames(X)
    rownames(ica_res$S) <- rownames(X)
    names(ica_res$stab) <- colnames(ica_res$M) <- colnames(ica_res$S) <- paste0("factor_", 1:ncol(ica_res$S))

    return(ica_res)
}

#' Stability ICA method
#' @import ica, BiocParallel
.stability_ica <- function(X, nc, resample, method, n_runs, BPPARAM, ...) {

    .ica_random <- function(i, nc, method, resample) {
        
        # Randomly initialises ICA
        set.seed(i)
        Rmat = matrix(rnorm(nc ** 2), nrow = nc, ncol = nc)

        if (resample) {
            X_bs <- X[, sample(ncol(X), replace = TRUE)]
        } else {
            X_bs <- X
        }
        
        # Get ICA loadings for given initialisation (and possibly bootstrap resample)
        S <- ica::ica(X_bs, nc=nc, method=method, center=FALSE, Rmat = Rmat, ...)$S
        colnames(S) <- paste0("seed_", i, "_", 1:nc)
        
        return(S)
    }

    S_all <- BiocParallel::bplapply(1:n_runs, .ica_random, BPPARAM = BPPARAM, nc=nc, method=method, resample=resample)
    S_all <- do.call(cbind, S_all)

    # Get correlations between factors and resulting clusters
    S_cor <- abs(cor(S_all))
    S_clust <- factor(cutree(hclust(as.dist(1 - S_cor)), k = nc))
    names(S_clust) <- colnames(S_all)

    stabilities <- c()
    centrotypes <- data.frame(matrix(nrow = nrow(S_all), ncol = nc, dimnames = list(rownames(S_all), 1:nc)))

    for (comp in 1:nc) {

        cluster_labels <- names(S_clust)[S_clust == comp]
        non_cluster_labels <- names(S_clust)[!names(S_clust) %in% cluster_labels]

        # Average intra-cluster similarity
        aics = mean(S_cor[cluster_labels, cluster_labels])

        # Average extra-cluster similarity
        aecs = mean(S_cor[cluster_labels, non_cluster_labels])

        stabilities <- c(stabilities, aics - aecs)

        which_is_centrotype <- which.max(apply(S_cor[cluster_labels, cluster_labels], 2, sum))
        centrotypes[[comp]] <- S_all[, cluster_labels[which_is_centrotype]]
    }

    idx <- order(stabilities, decreasing = TRUE)
    centrotypes <- centrotypes[, idx]
    stabilities <- stabilities[idx]

    return(list(stab = stabilities, S = centrotypes))
}

#' @import moments
.reorient_factors <- function(S) {
    skew <- ifelse(apply(S, 2, moments::skewness) >= 0, 1, -1)
    for (i in 1:ncol(S)) {
        S[,i] <- S[,i] * skew[i]
    }
    return(S)
}

#' X = M * S
#' X / S = M
#' X * inv(S) = M
#' @import MASS
.project_ica <- function(newdata, S) {
    M <- t(newdata) %*% t(MASS::ginv(S))
    colnames(M) <- colnames(S)

    return(M)
}

estimate_components <- function(X, min_components=10, max_components=60, by=2,
                                min_mean_stability = 0.85, n_runs = 30,
                                resample = FALSE,
                                center_X=TRUE, scale_X=FALSE, verbose = TRUE, 
                                BPPARAM = SerialParam(), ...) {
    if (inherits(X, "SummarizedExperiment")) {
        X <- assay(X, "normal")
    }

    stabilities <- data.frame()
    
    if (verbose) tpb <- txtProgressBar(min = min_components, max = max_components, initial = min_components, style = 3)

    for (nc in seq(from = min_components, to = max_components, by = by)) {

        ica_res <- run_ica(X, nc=nc, center_X=center_X, scale_X=scale_X, use_stability=TRUE, resample=resample, BPPARAM=BPPARAM, method="fast", n_runs=n_runs, ...)

        stabilities <- rbind(stabilities, data.frame(
            nc = nc,
            comp = names(ica_res$stab),
            stability = ica_res$stab
        ))
        
        if (verbose) setTxtProgressBar(tpb, nc)
    }
    
    if (verbose) close(tpb)

    mean_stabilities <- aggregate(stabilities, list(stabilities$nc), mean)
    select_nc <- min(mean_stabilities$nc[mean_stabilities$stability >= min_mean_stability])

    return(list(stabilities = stabilities, select_nc = select_nc))
}
