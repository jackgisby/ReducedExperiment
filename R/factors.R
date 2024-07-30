#' Apply dimensionality reduction using ICA
#'
#' Performs independent component analysis (ICA).
#' Calls \link{ReducedExperiment}[run_ica] to perform the analysis.
#'
#' @param X Either a \link{SummarizedExperiment}[SummarizedExperiment] object
#' or a matrix containing data to be subject to ICA. `X` should have rows as
#' features and columns as samples.
#'
#' @param nc The number of components to be identified. See
#' \link{ReducedExperiment}[estimate_stability] for a method to estimate the
#' optimal number of components.
#'
#' @param center_X If TRUE, X is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to ICA. Generally recommended.
#'
#' @param scale_X If TRUE, X is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before ICA.
#'
#' @param assay_name If `X` is a
#' \link{SummarizedExperiment}[SummarizedExperiment], then this should be the
#' name of the assay to be subject to ICA.
#'
#' @param ... Additional arguments to be passed to
#' \link{ReducedExperiment}[run_ica].
#'
#' @export
estimate_factors <- function(
        X,
        nc,
        center_X = TRUE,
        scale_X = FALSE,
        assay_name = "normal",
        ...) {
    if (!inherits(X, "SummarizedExperiment")) {
        X <- SummarizedExperiment(assays = list("normal" = X))
    }

    if (assay_name != "normal") {
        assay(X, "normal") <- assay(X, assay_name)
    }

    if ("transformed" %in% assayNames(X)) {
        warning("Overwriting 'transformed' assay slot in X")
    }

    assay(X, "transformed") <- t(scale(t(assay(X, "normal")),
        center = center_X,
        scale = scale_X
    ))

    # print(assay(X, "transformed"))
    if (center_X) center_X <- attr(assay(X, "transformed"), "scaled:center")
    if (scale_X) scale_X <- attr(assay(X, "transformed"), "scaled:scale")

    ica_res <- run_ica(assay(X, "transformed"),
        nc = nc,
        center_X = FALSE, scale_X = FALSE, ...
    )

    return(.se_to_fe(
        X,
        reduced = ica_res$M,
        loadings = ica_res$S,
        stability = ica_res$stab,
        center_X = center_X,
        scale_X = scale_X
    ))
}

#' Creates a FactorisedExperiment from a SummarizedExperiment
#'
#' Helper function for transforming a
#' \link{ReducedExperiment}[FactorisedExperiment] into a
#' \link{SummarizedExperiment}[SummarizedExperiment]
#'
#' @param se A \link{SummarizedExperiment}[SummarizedExperiment] object.
#'
#' @param reduced Data to be passed to the `reduced` slot.
#'
#' @param loadings Data to be passed to the `loadings` slot.
#'
#' @param stability Data to be passed to the `stability` slot.
#'
#' @param center_X Data to be passed to the `center_X` slot.
#'
#' @param scale_X Data to be passed to the `scale_X` slot.
#'
#' @noRd
#' @keywords internal
.se_to_fe <- function(se, reduced, loadings, stability, center_X, scale_X) {
    return(FactorisedExperiment(
        loadings = loadings, stability = stability,
        center = center_X, scale = scale_X, reduced = reduced,
        assays = assays(se), rowData = rowData(se),
        colData = colData(se), metadata = S4Vectors::metadata(se)
    ))
}

#' Run ICA for a data matrix
#'
#' Runs ICA through \link{ica}[ica]. X is passed directly to \link{ica}[ica],
#' with rows as features and samples as columns. Generates a source signal
#' matrix (loadings) with rows as features and columns as factors and a mixture
#' matrix (reduced data) with samples as rows and columns as factors.
#'
#' @param X A matrix with features as rows and columns as samples.
#'
#' @param nc The number of components to be identified. See
#' \link{ReducedExperiment}[estimate_stability] for a method to estimate the
#' optimal number of components.
#'
#' @param use_stability Whether to use a stability-based approach to estimate
#' factors. See `details` for further information
#'
#' @param resample If TRUE, a boostrap approach is used to estimating factors.
#' Else, random initialisation of ICA is employed. Ignored if `use_stability`
#' is FALSE. See `details` for further information.
#'
#' @param method The ICA method to use. Passed to \link{ica}[ica], the options
#' are "fast", "imax" or "jade".
#'
#' @param stability_threshold A stability threshold for pruning factors. Factors
#' with a stability below this threshold will be removed. If used, the threshold
#' can lead to fewer factors being returned than that specified by `nc`.
#'
#' @param center_X If TRUE, X is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to ICA. Generally recommended.
#'
#' @param scale_X If TRUE, X is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before ICA.
#'
#' @param reorient_skewed If TRUE, factors are reorientated to ensure that the
#' loadings of each factor (i.e., the source signal matrix) have positive skew.
#' Helps ensure that the most influential features for each factor are
#' positively associated with it.
#'
#' @param scale_components If TRUE, the loadings are standardised (to have a
#' mean of 0 and standard deviation of 1).
#'
#' @param scale_reduced If TRUE, the reduced data (mixture matrix) are
#' standardised (to have a mean of 0 and standard deviation of 1).
#'
#' @param n_runs The number of times to run ICA. Ignored if `use_stability` is
#' FALSE. See `details` for further information.
#'
#' @param BPPARAM A class containing parameters for parallel evaluation. Uses
#' \link[BiocParallel]{SerialParam} by default, running only a single
#' ICA computation at a time. Ignored if `use_stability`
#' is FALSE. See `details` for further information.
#'
#' @param ... Additional arguments to be passed to
#' \link{ica}[ica].
#'
#' @details
#' Function performs ICA for a data matrix. If `use_stability` is TRUE, then
#' ICA is performed multiple times with either: i) random initialisation
#' (default); or ii) bootstrap resampling of the data (if `resample` is TRUE).
#'
#' The stability-based ICA algorithm is similar to the the ICASSO approach
#' (https://www.cs.helsinki.fi/u/ahyvarin/papers/Himberg03.pd) that is
#' implemented in the stabilized-ica Python package
#' (https://github.com/ncaptier/stabilized-ica/tree/master).
#'
#' Results from this function should be broadly similar to those generated by
#' stabilized-ica, although they will not be identical. Notable differences
#' include:
#' \describe{
#'  \item{"ICA algorithm"}{Differences in the underlying implementation of
#'  ICA could contribute to changes in the final results.}
#'  \item{"Stability threshold"}{The `stability_threshold` argument, if
#'  specified, removes unstable components. Such a threshold is not
#'  used by stabilized-ica.}
#'  \item{"Mixture matrix recovery"}{ICA is generally formulated as
#'  `X = MS`, where `X` is the input data, `M` is the mixture matrix
#'  (reduced data) and `S` is the source signal matrix (feature loadings).
#'  The stabilised ICA approach first calculates a source signal matrix
#'  before recovering the mixture matrix. To do this, other implementations,
#'  including that of the stabilized-ica package, multiply `X` by the
#'  pseudo-inverse of `S`. Such an operation is implemented in the `ginv`
#'  function of the `MASS` R package. In the development of ReducedExperiment,
#'  we noticed that taking the inverse of `S` often failed, particularly when
#'  there were correlated factors. For this reason, we instead formulate the
#'  mixture matrix as `M = XS`. After standardisation of `M`, both approaches
#'  return near-identical results, given that the matrix inverse was
#'  successfully calculated.}
#' }
#'
#' In short, the stability-based algorithm consists of:
#' i) Running ICA multiple times with either random initialisation or bootstrap
#' resampling of the input data.
#' ii) Clustering the resulting factors across all runs based on the
#' signature matrix.
#' iii) Calculating intra- (aics) and extra- (aecs) cluster
#' stability, and defining the final cluster stability as `aics - aecs`.
#' iv) Calculating the cluster centrotype as the factor with the highest
#' intra-cluster stability.
#' v) Optionally removing factors below a specified stability threshold
#' (`stability_threshold`).
#'
#' @import ica
#' @export
run_ica <- function(
        X, nc, use_stability = FALSE, resample = FALSE,
        method = "fast", stability_threshold = NULL,
        center_X = TRUE, scale_X = FALSE,
        reorient_skewed = TRUE,
        scale_components = TRUE, scale_reduced = TRUE,
        n_runs = 30,
        BPPARAM = BiocParallel::SerialParam(RNGseed = 1),
        ...) {

    if (center_X | scale_X) {
        X <- t(scale(t(X), center = center_X, scale = scale_X))
    }

    if (use_stability) {
        ica_res <- .stability_ica(
            X,
            nc = nc,
            resample = resample,
            method = method,
            stability_threshold = stability_threshold,
            n_runs = n_runs,
            BPPARAM = BPPARAM,
            ...
        )
    } else {
        if (resample) {
            stop("Cannot use resampling approach when `use_stability` is FALSE")
        }

        if (!is.null(stability_threshold)) {
            stop("Cannot apply `stability_threshold` when `use_stability` is
                 FALSE")
        }

        ica_res <- list(S = ica::ica(
            X,
            nc = nc,
            method = method,
            center = FALSE,
            ...
        )$S)
    }

    # Reorient and scale factors before recalculating M
    if (reorient_skewed) ica_res$S <- .reorient_factors(ica_res$S)
    if (scale_components) ica_res$S <- scale(ica_res$S)
    ica_res$M <- .project_ica(X, ica_res$S)
    if (scale_reduced) ica_res$M <- scale(ica_res$M)

    # Add factors / sample names
    rownames(ica_res$M) <- colnames(X)
    rownames(ica_res$S) <- rownames(X)
    colnames(ica_res$M) <- colnames(ica_res$S) <-
        paste0("factor_", seq_len(ncol(ica_res$S)))

    if (use_stability) names(ica_res$stab) <- colnames(ica_res$S)

    return(ica_res)
}

#' Stability ICA method
#'
#' Function for running stabilised ICA. See \link{ReducedExperiment}[run_ica].
#'
#' @import ica
#' @import BiocParallel
#'
#' @noRd
#' @keywords internal
.stability_ica <- function(
        X,
        nc,
        resample,
        method,
        n_runs,
        BPPARAM,
        stability_threshold,
        BPOPTIONS = bpoptions(),
        return_centrotypes = TRUE,
        ...) {
    .ica_random <- function(i, nc, method, resample) {

        # Randomly initialises ICA
        Rmat <- matrix(stats::rnorm(nc**2), nrow = nc, ncol = nc)

        if (resample) {
            X_bs <- X[, sample(ncol(X), replace = TRUE)]
        } else {
            X_bs <- X
        }

        # Get ICA loadings for given initialisation
        # (and possibly bootstrap resample)
        S <- ica::ica(
            X_bs,
            nc = nc,
            method = method,
            center = FALSE,
            Rmat = Rmat,
            ...
        )$S

        colnames(S) <- paste0("iteration_", i, "_", seq_len(ncol(S)))
        if (ncol(S) != nc) warning("ICA did not return expected number of
                                   factors, potentially indicating a rank
                                   deficiency in the input")

        return(S)
    }

    S_all <- BiocParallel::bplapply(
        seq_len(n_runs),
        .ica_random,
        BPPARAM = BPPARAM,
        nc = nc,
        method = method,
        resample = resample
    )
    S_all <- do.call(cbind, S_all)

    # Get correlations between factors and resulting clusters
    S_cor <- abs(stats::cor(S_all))
    S_clust <- factor(stats::cutree(stats::hclust(stats::as.dist(1 - S_cor)),
        k = nc
    ))
    names(S_clust) <- colnames(S_all)

    stabilities <- c()
    centrotypes <- data.frame(
        matrix(
            nrow = nrow(S_all),
            ncol = nc,
            dimnames = list(rownames(S_all), seq_len(nc))
        )
    )

    for (comp in seq_len(nc)) {
        cluster_labels <- names(S_clust)[S_clust == comp]
        non_cluster_labels <-
            names(S_clust)[!names(S_clust) %in% cluster_labels]

        # Average intra-cluster similarity
        aics <- mean(S_cor[cluster_labels, cluster_labels])

        # Average extra-cluster similarity
        aecs <- mean(S_cor[cluster_labels, non_cluster_labels])

        stabilities <- c(stabilities, aics - aecs)

        which_is_centrotype <-
            which.max(apply(S_cor[cluster_labels, cluster_labels], 2, sum))
        centrotypes[[comp]] <- S_all[, cluster_labels[which_is_centrotype]]
    }

    if (!return_centrotypes) {
        return(list(
            stab = stabilities,
            S_all = S_all,
            S_clust = S_clust,
            S_cor = S_cor
        ))
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

#' Ensure factors have positive skew
#'
#' Reorientates factors based on their skew. Generally ensures that the most
#' aligned genes are positively associated with the factor.
#'
#' @import moments
#'
#' @noRd
#' @keywords internal
.reorient_factors <- function(S) {
    skew <- ifelse(apply(S, 2, moments::skewness) >= 0, 1, -1)
    for (i in seq_len(ncol(S))) {
        S[, i] <- S[, i] * skew[i]
    }
    return(S)
}

#' Calculate the mixture matrix
#'
#' Calculates the mixture matrix (reduced data, sample-level) from the source
#' signal matrix (loadings, gene-level) and original data
#' (samples vs. features).
#'
#' Originally used MASS:ginv(S) to do the following calculation:
#'
#' X = M * S
#' X / S = M
#' X * inv(S) = M
#'
#' Now instead does the following:
#' M = X * S
#'
#' @noRd
#' @keywords internal
.project_ica <- function(newdata, S) {
    # M <- t(newdata) %*% t(MASS::ginv(S))
    M <- t(newdata) %*% S

    colnames(M) <- colnames(S)

    return(M)
}

#' Estimate stability of factors as a function of the number of components
#'
#' Estimates the stability of factors over a range of component numbers to
#' aid in the identification of the optimal factor number.
#'
#' @param X Either a \link{SummarizedExperiment}[SummarizedExperiment] object
#' or a matrix containing data to be subject to ICA. `X` should have rows as
#' features and columns as samples.
#'
#' @param min_components The minimum number of components to estimate the
#' stability for.
#'
#' @param max_components The maximum number of components to estimate the
#' stability for.
#'
#' @param by The number by which to increment the numbers of components
#' tested.
#'
#' @param n_runs The number of times to run ICA. Ignored if `use_stability` is
#' FALSE. See `details` for further information.
#'
#' @param resample If TRUE, a boostrap approach is used to estimating factors.
#' Else, random initialisation of ICA is employed. Ignored if `use_stability`
#' is FALSE. See `details` for further information.
#'
#' @param mean_stability_threshold The function will return the minimal
#' number of components that exceed this stability threshold.
#'
#' @param center_X If TRUE, X is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to ICA. Generally recommended.
#'
#' @param scale_X If TRUE, X is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before ICA.
#'
#' @param assay_name If `X` is a
#' \link{SummarizedExperiment}[SummarizedExperiment], then this should be the
#' name of the assay to be subject to ICA.
#'
#' @param BPPARAM A class containing parameters for parallel evaluation. Uses
#' \link[BiocParallel]{SerialParam} by default, running only a single
#' ICA computation at a time. Ignored if `use_stability`
#' is FALSE. See `details` for further information.
#'
#' @param verbose If TRUE, shows a progress bar that updates for each
#' number of components tested. Note that the time taken may not be linear,
#' because the time taken to run ICA generally increases with the number
#' of components.
#'
#' @param ... Additional arguments to be passed to
#' \link{ReducedExperiment}[run_ica].
#'
#' @details
#' Runs the stability-based ICA algorithm
#' (see \link[ReducedExperiment]{run_ica}) for a range of component numbers.
#' Estimates stability for each one, allowing for selection of the optimal
#' number of components to be used for ICA. The results of this function
#' can be plotted by \link[ReducedExperiment]{plot_stability}.
#'
#' This algorithm is similar to the Most Stable Transcriptome
#' Dimension (MSTD) approach
#' (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4112-9).
#'
#' The function automatically selects a number of components based on
#' `mean_stability_threshold`. Ideally, this should not be trusted blindly,
#' but instead a choice can be made based on visualisation of the
#' stabilities, as created by \link[ReducedExperiment]{plot_stability}. The
#' MSTD paper provides additional context and advice on choosing the number
#' of components based on these data.
#'
#' @returns Returns a list containing: i) a list of stabilities as a function
#' of the number of components; and ii) the selected number of components based
#' on the `mean_stability_threshold`.
#'
#' @export
estimate_stability <- function(
    X,
    min_components = 10,
    max_components = 60,
    by = 2,
    n_runs = 30,
    resample = FALSE,
    mean_stability_threshold = NULL,
    center_X = TRUE,
    scale_X = FALSE,
    assay_name = "normal",
    BPPARAM = BiocParallel::SerialParam(RNGseed = 1),
    verbose = TRUE,
    ...
) {
    if (inherits(X, "SummarizedExperiment")) {
        X <- assay(X, "normal")
    }

    if (assay_name != "normal") {
        assay(X, "normal") <- assay(X, assay_name)
    }

    if (dim(X)[2] < max_components) stop("Number of samples must be greater
                                         than max_components")

    stabilities <- data.frame()

    if (verbose) {
        tpb <- utils::txtProgressBar(
            min = min_components,
            max = max_components,
            initial = min_components,
            style = 3
        )
    }

    for (nc in seq(from = min_components, to = max_components, by = by)) {
        ica_res <- run_ica(
            X,
            nc = nc,
            center_X = center_X,
            scale_X = scale_X,
            use_stability = TRUE,
            resample = resample,
            BPPARAM = BPPARAM,
            method = "fast",
            n_runs = n_runs,
            ...
        )

        stabilities <- rbind(stabilities, data.frame(
            nc = nc,
            component_name = names(ica_res$stab),
            component_number = as.numeric(gsub(
                "factor_", "",
                names(ica_res$stab)
            )),
            stability = ica_res$stab
        ))

        if (verbose) utils::setTxtProgressBar(tpb, nc)
    }

    if (verbose) close(tpb)

    select_nc <- NULL

    if (!is.null(mean_stability_threshold)) {
        mean_stabilities <-
            stats::aggregate(stabilities$stability, list(stabilities$nc), mean)
        colnames(mean_stabilities) <- c("nc", "stability")

        if (any(mean_stabilities$stability >= mean_stability_threshold)) {
            select_nc <- max(mean_stabilities$nc[mean_stabilities$stability >=
                mean_stability_threshold])
        }
    }

    return(list("stability" = stabilities, "selected_nc" = select_nc))
}

#' Plot component stability as a function of the number of components
#'
#' Plots the results of \link[ReducedExperiment]{estimate_stability}.
#'
#' @param stability The results of \link[ReducedExperiment]{estimate_stability}.
#'
#' @param plot_path The path at which the plot will be saved
#'
#' @param stability_threshold Plots a stability threshold, below which
#' components can be pruned by \link[ReducedExperiment]{run_ica}.
#'
#' @param mean_stability_threshold Plots a stability threshold, which is used
#' by \link[ReducedExperiment]{estimate_stability} to estimate the optimal
#' number of components.
#'
#' @param height The height of the plot, to be passed to \link[ggplot2]{ggsave}.
#'
#' @param width The width of the plot, to be passed to \link[ggplot2]{ggsave}.
#'
#' @param ... Additional arguments to be passed to \link[ggplot2]{ggsave}.
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
plot_stability <- function(
        stability, plot_path = NULL,
        stability_threshold = NULL, mean_stability_threshold = NULL,
        height = 4, width = 10, ...) {
    if (is.list(stability)) stability <- stability[["stability"]]

    stab_plot <- ggplot(stability, aes(
            !!sym("component_number"),
            !!sym("stability"),
            group = !!sym("nc")
        )) +
        geom_line() +
        ylim(c(0, 1)) +
        ylab("Component stability") +
        xlab("Component number")

    if (!is.null(stability_threshold)) {
        stab_plot <- stab_plot +
            geom_hline(yintercept = stability_threshold)
    }

    stabilities_agg <- stats::aggregate(
        stability[, c("nc", "stability")],
        list(stability$nc),
        mean
    )

    mean_stab_plot <- ggplot(stabilities_agg,
                             aes(!!sym("nc"), !!sym("stability"), group = 1))  +
        geom_line() +
        ylim(c(0, 1)) +
        ylab("Mean component stability") +
        xlab("Number of components")

    if (!is.null(mean_stability_threshold)) {
        mean_stab_plot <- mean_stab_plot +
            geom_hline(yintercept = mean_stability_threshold)
    }

    combined_plot <- stab_plot + mean_stab_plot

    if (!is.null(plot_path)) {
        ggsave(plot_path, combined_plot,
            height = height, width = width, ...
        )
    }

    return(list(
        "combined_plot" = combined_plot,
        "stability_plot" = stab_plot,
        "mean_plot" = mean_stab_plot
    ))
}
