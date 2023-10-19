#' Apply dimensionality reduction
#' @export
reduce_data <- function(
        X, method = c("ICA", "WGCNA"),
        phenoData=Biobase::annotatedDataFrameFrom(X, byrow=FALSE),
        featureData=Biobase::annotatedDataFrameFrom(X, byrow=TRUE),
        center_X=(method == "ICA"), scale_X=FALSE,
        experimentData=MIAME(), annotation=character(),
        protocolData=annotatedDataFrameFrom(X, byrow=FALSE),
        ...)
{
    if ("ExpressionSet" %in% class(X)) {
        message("Using data from ExpressionSet and ignoring other arguments.")

        phenoData <- phenoData(X)
        featureData <- featureData(X)
        experimentData <- experimentData(X)
        annotation <- annotation(X)
        protocolData <- protocolData(X)
        X <- exprs(X)
    }

    X_tr <- scale(t(X), center=center_X, scale=scale_X)

    if (method == "ICA") {
        ica_res <- run_ica(X_tr, center_X=FALSE, ...)
        reduced_set <- FactorSet(exprs=t(X_tr), reduced=ica_res$M, S=ica_res$S,
                                 varexp=ica_res$vafs, phenoData=phenoData,
                                 featureData=featureData, annotation=annotation,
                                 experimentData=experimentData,
                                 protocolData=protocolData)

    } else if (method == "WGCNA") {
        wgcna_res <- run_wgcna(X, ...)
        reduced_set <- ModuleSet(X, wgcna_res$E, wgcna_res$assignments,
                                 phenoData=phenoData, featureData=featureData,
                                 annotation=annotation,
                                 experimentData=experimentData,
                                 protocolData=protocolData)
    } else {
        stop(paste0("Method ", method, " not recognised."))
    }

    # TODO: Store these attributes somewhere so new data can be projected
    # if (!is.null(attr(x, "scaled:scale")) | !is.null(attr(x, "scaled:center")))
    #     reduced_set

    return(reduced_set)
}

#' Project data
#' @export
project_data <- function(reduced_set, newdata) {

}

#' Run ICA
#' @export
run_ica <- function(X, nc, method="imax", maxit=500, tol=1e-6, seed=1,
                    center_X=TRUE, scale_X=FALSE, reorient_skewed=TRUE, ...) {
    set.seed(seed)

    if (center_X | scale_X)
        {X <- scale(X, center=center_X, scale=scale_X)}

    ica_res <- ica::ica(t(X), nc = nc, maxit = maxit, tol = tol,
                        method = method, ... )

    # Add factors / sample names
    rownames(ica_res$M) <- rownames(X)
    colnames(ica_res$M) <- colnames(ica_res$S) <- paste0("factor_",
                                                         1:ncol(ica_res$S))

    if (reorient_skewed) {ica_res$S <- .reorient_skewed(ica_res$S)}

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
    M <- newdata %*% t(MASS::ginv(S))
    colnames(M) <- colnames(S)

    return(M)
}

#' Run WGCNA
#' @export
run_wgcna <- function() {
    stop("Not implemented yet")
}
