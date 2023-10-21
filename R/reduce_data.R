#' Apply dimensionality reduction using ICA
#' @export
run_ica <- function(X, nc, method="imax", center_X=TRUE, scale_X=FALSE,
                    reorient_skewed=TRUE, seed=1, ...)
{
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list(normal = X))
    }

    assay(X, "transformed") <- t(scale(t(assay(X, "normal")), center=center_X, scale=scale_X))

    ica_res <- compute_ica(assay(X, "transformed"), nc=nc, method=method,
                           center_X=FALSE,  reorient_skewed=reorient_skewed,
                           seed=seed, ...)

    return(.se_to_fe(X, reduced=ica_res$M, loadings=ica_res$S, varexp=ica_res$vafs))
}

.se_to_fe <- function(se, reduced, loadings, varexp) {
    return(FactorisedExperiment(loadings=loadings, varexp=varexp, reduced=reduced,
                                assays=assays(se), rowData=rowData(se),
                                colData=colData(se), metadata=metadata(se)))
}

.se_to_me <- function(se, reduced, assignments) {
    return(FactorisedExperiment(assignments=assignments, reduced=reduced,
                                assays=assays(se), rowData=rowData(se),
                                colData=colData(se), metadata=metadata(se)))
}

#' Run ICA for a data matrix
#' @export
compute_ica <- function(X, nc, method="imax", center_X=TRUE, scale_X=FALSE,
                        reorient_skewed=TRUE, seed=1, ...) {
    set.seed(seed)

    if (center_X | scale_X)
        {X <- scale(X, center=center_X, scale=scale_X)}

    ica_res <- ica::ica(X, nc=nc, method = method, ... )

    # Reorient factors and recalculate M
    if (reorient_skewed) {ica_res$S <- .reorient_skewed(ica_res$S)}
    ica_res$M <- .project_ica(X, ica_res$S)

    # Add factors / sample names
    rownames(ica_res$M) <- colnames(X)
    rownames(ica_res$S) <- rownames(X)
    names(ica_res$vafs) <- colnames(ica_res$M) <- colnames(ica_res$S) <- paste0("factor_", 1:ncol(ica_res$S))

    return(ica_res)
}

.reorient_skewed <- function(S) {
    skew <- ifelse(apply(S, 2, moments::skewness) >= 0, 1, -1)
    for (i in 1:ncol(S)) {
        S[,i] <- S[,i] * skew[i]
    }
    return(S)
}

#' X = M * S
#' X / S = M
#' X * inv(S) = M
.project_ica <- function(newdata, S) {
    M <- t(newdata) %*% t(MASS::ginv(S))
    colnames(M) <- colnames(S)

    return(M)
}

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
