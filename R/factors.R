#' Apply dimensionality reduction using ICA
#' @export
estimate_factors <- function(X, nc, center_X=TRUE, scale_X=FALSE, assay_name="normal", ...)
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

    ica_res <- run_ica(assay(X, "transformed"), nc=nc,
                       center_X=FALSE, scale_X=FALSE, ...)

    return(.se_to_fe(X, reduced=ica_res$M, loadings=ica_res$S, stability=ica_res$stab, center_X=center_X, scale_X=scale_X))
}

.se_to_fe <- function(se, reduced, loadings, stability, center_X, scale_X) {
    return(FactorisedExperiment(loadings=loadings, stability=stability,
                                center=center_X, scale=scale_X, reduced=reduced,
                                assays=assays(se), rowData=rowData(se),
                                colData=colData(se), metadata=metadata(se)))
}

#' Run ICA for a data matrix
#' @import ica
#' @export
run_ica <- function(X, nc, use_stability=FALSE, resample=FALSE,
                    method="fast", stability_threshold=NULL,
                    center_X=TRUE, scale_X=FALSE,
                    reorient_skewed=TRUE, seed=1,
                    scale_components=TRUE, scale_reduced=TRUE,
                    n_runs=30,
                    BPPARAM = BiocParallel::SerialParam(),
                    ...) {
    set.seed(seed)

    if (center_X | scale_X)
        {X <- t(scale(t(X), center=center_X, scale=scale_X))}

    if (use_stability) {
        ica_res <- .stability_ica(X, nc=nc, resample=resample, method=method, stability_threshold=stability_threshold, n_runs=n_runs, BPPARAM=BPPARAM, ...)
    } else {
        if (resample) stop("Cannot use resampling approach when `use_stability` is FALSE")
        if (!is.null(stability_threshold)) stop("Cannot apply `stability_threshold` when `use_stability` is FALSE")

        ica_res <- list(S = ica::ica(X, nc=nc, method=method, center=FALSE, ...)$S)
    }

    # Reorient and scale factors before recalculating M
    if (reorient_skewed) ica_res$S <- .reorient_factors(ica_res$S)
    if (scale_components) ica_res$S <- scale(ica_res$S)
    ica_res$M <- .project_ica(X, ica_res$S)
    if (scale_reduced) ica_res$M <- scale(ica_res$M)

    # Add factors / sample names
    rownames(ica_res$M) <- colnames(X)
    rownames(ica_res$S) <- rownames(X)
    colnames(ica_res$M) <- colnames(ica_res$S) <- paste0("factor_", 1:ncol(ica_res$S))
    if (use_stability) names(ica_res$stab) <- colnames(ica_res$S)

    return(ica_res)
}

#' Stability ICA method
#' @import ica
#' @import BiocParallel
.stability_ica <- function(X, nc, resample, method, n_runs, BPPARAM, stability_threshold,
                           BPOPTIONS = bpoptions(), return_centrotypes = TRUE, ...) {

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
        colnames(S) <- paste0("seed_", i, "_", 1:ncol(S))

        if (ncol(S) != nc) warning("ICA did not return expected number of factors, potentially indicating a rank deficiency in the input")

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

    if (!return_centrotypes) {
        return(list(stab = stabilities, S_all = S_all, S_clust = S_clust, S_cor = S_cor))
    }

    stability_order <- order(stabilities, decreasing = TRUE)
    centrotypes <- centrotypes[, stability_order]
    stabilities <- stabilities[stability_order]

    if (!is.null(stability_threshold)) {
        above_stability_threshold <- which(stabilities > stability_threshold)
        centrotypes <- centrotypes[, above_stability_threshold]
        stabilities <- stabilities[above_stability_threshold]
    }

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
.project_ica <- function(newdata, S) {

    # M <- t(newdata) %*% t(MASS::ginv(S))
    M <- t(newdata) %*% S

    colnames(M) <- colnames(S)

    return(M)
}

estimate_stability <- function(X, min_components=10, max_components=60,
                                    by=2, n_runs = 30, resample = FALSE,
                                    mean_stability_threshold = NULL,
                                    center_X=TRUE, scale_X=FALSE,
                                    BPPARAM = BiocParallel::SerialParam(),
                                    verbose = TRUE, ...) {
    if (inherits(X, "SummarizedExperiment")) {
        X <- assay(X, "normal")
    }

    stabilities <- data.frame()

    if (verbose) tpb <- txtProgressBar(min = min_components, max = max_components, initial = min_components, style = 3)

    for (nc in seq(from = min_components, to = max_components, by = by)) {

        ica_res <- run_ica(X, nc=nc, center_X=center_X, scale_X=scale_X, use_stability=TRUE, resample=resample, BPPARAM=BPPARAM, method="fast", n_runs=n_runs, ...)

        stabilities <- rbind(stabilities, data.frame(
            nc = nc,
            component_name = names(ica_res$stab),
            component_number = as.numeric(gsub("factor_", "", names(ica_res$stab))),
            stability = ica_res$stab
        ))

        if (verbose) setTxtProgressBar(tpb, nc)
    }

    if (verbose) close(tpb)

    select_nc <- NULL

    if (!is.null(mean_stability_threshold)) {
        mean_stabilities <- aggregate(stabilities$stability, list(stabilities$nc), mean)
        colnames(mean_stabilities) <- c("nc", "stability")

        if (any(mean_stabilities$stability >= mean_stability_threshold)) {
            select_nc <- max(mean_stabilities$nc[mean_stabilities$stability >= mean_stability_threshold])
        }
    }

    return(list("stability" = stabilities, "selected_nc" = select_nc))
}

#' Plot component stability as a function of the number of components
#'
#' @references MSTD
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
plot_stability <- function(stability, plot_path,
                           stability_threshold=NULL, mean_stability_threshold=NULL,
                           height = 4, width = 10, ...) {

    if (is.list(stability)) stability <- stability[["stability"]]

    stab_plot <- ggplot(stability, aes(component_number, stability, group = nc)) +
        geom_line() +
        ylim(c(0,1)) +
        ylab("Component stability") +
        xlab("Component number")

    if (!is.null(stability_threshold)) {
        stab_plot <- stab_plot +
            geom_hline(yintercept = stability_threshold)
    }

    mean_stab_plot <- ggplot(aggregate(stability, list(stability$nc), mean), aes(nc, stability, group = 1)) +
        geom_line() +
        ylim(c(0,1)) +
        ylab("Mean component stability") +
        xlab("Number of components")

    if (!is.null(mean_stability_threshold)) {
        mean_stab_plot <- mean_stab_plot +
            geom_hline(yintercept = mean_stability_threshold)
    }

    combined_plot <- stab_plot + mean_stab_plot

    if (!is.null(plot_path)) ggsave(plot_path, combined_plot, height=height, width=width, ...)

    return(list(
        "combined_plot" = combined_plot,
        "stability_plot" = stab_plot,
        "mean_plot" = mean_stab_plot
    ))
}
