#' Perform dimensionality reduction using Independent Component Analysis
#'
#' Performs independent component analysis (ICA) and packages both the input
#' data and subsequent results into a
#' \link[ReducedExperiment]{FactorisedExperiment} container. Calls
#' \link[ReducedExperiment]{runICA} to perform the analysis; see its
#' documentation page for more information on the ICA method, parameters
#' and outputs.
#'
#' @param X Either a \link[SummarizedExperiment]{SummarizedExperiment} object
#' or a matrix containing data to be subject to ICA. `X` should have rows as
#' features and columns as samples.
#'
#' @param nc The number of components to be identified. See
#' \link[ReducedExperiment]{estimateStability} for a method to estimate the
#' optimal number of components.
#'
#' @param center_X If `TRUE`, X is centered (i.e., features / rows are transformed
#' to have a mean of 0) prior to ICA. Generally recommended.
#'
#' @param scale_X If `TRUE`, X is scaled (i.e., features / rows are transformed
#' to have a standard deviation of 1) before ICA.
#'
#' @param assay_name If `X` is a
#' \link[SummarizedExperiment]{SummarizedExperiment}, then this should be the
#' name of the assay to be subject to ICA.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{runICA}.
#'
#' @returns A \link[ReducedExperiment]{FactorisedExperiment} is returned
#' containing the input data (i.e., the original data matrix in addition to
#' other slots if a \link[SummarizedExperiment]{SummarizedExperiment} was used
#' as input). Additionally contains the results of factor analysis, stored in
#' the `reduced` and `loadings` slots. The `center_X`, `scale_X` and
#' `stability` slots may also be filled depending on the arguments given
#' to `estimateFactors`.
#'
#' @seealso [ReducedExperiment::runICA()], [ica::ica()]
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 100 rows (features)
#' # and 20 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(100, 20, "feature", "obs")
#'
#' # Estimate 5 factors based on the data matrix
#' set.seed(1)
#' fe_1 <- estimateFactors(X, nc = 5)
#' fe_1
#'
#' # Convert the data matrix to a SummarizedExperiment, then estimate 5 factors
#' se <- SummarizedExperiment(assays = list("normal" = X))
#' set.seed(1)
#' fe_2 <- estimateFactors(se, nc = 5)
#' fe_2
#'
#' @export
estimateFactors <- function(
    X,
    nc,
    center_X = TRUE,
    scale_X = FALSE,
    assay_name = "normal",
    ...
) {
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

    if (center_X) center_X <- attr(assay(X, "transformed"), "scaled:center")
    if (scale_X) scale_X <- attr(assay(X, "transformed"), "scaled:scale")

    ica_res <- runICA(assay(X, "transformed"),
        nc = nc,
        center_X = FALSE, scale_X = FALSE, ...
    )

    return(.seToFe(
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
#' \link[ReducedExperiment]{FactorisedExperiment} into a
#' \link[SummarizedExperiment]{SummarizedExperiment}
#'
#' @param se A \link[SummarizedExperiment]{SummarizedExperiment} object.
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
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.seToFe <- function(se, reduced, loadings, stability, center_X, scale_X) {
    return(FactorisedExperiment(
        loadings = loadings, stability = stability,
        center = center_X, scale = scale_X, reduced = reduced,
        assays = assays(se), rowData = rowData(se),
        colData = colData(se), metadata = S4Vectors::metadata(se)
    ))
}

#' Run standard or stabilised Independent Component Analysis
#'
#' Runs ICA through \link[ica]{ica}. If `use_stability` is FALSE, then `X` is
#' passed directly to \link[ica]{ica} and a standard ICA analysis is performed.
#' If `use_stability` is `TRUE`, then the stabilised ICA procedure is carried
#' out (see `details`).
#'
#' @param use_stability Whether to use a stability-based approach to estimate
#' factors. See `details` for further information.
#'
#' @param resample If `TRUE`, a boostrap approach is used to estimate factors
#' and quantify stability. Else, random initialisation of ICA is employed.
#' Ignored if `use_stability` is `FALSE`.
#'
#' @param method The ICA method to use. Passed to \link[ica]{ica}, the options
#' are "fast", "imax" or "jade".
#'
#' @param stability_threshold A stability threshold for pruning factors. Factors
#' with a stability below this threshold will be removed. If used, the threshold
#' can lead to fewer factors being returned than that specified by `nc`.
#'
#' @param reorient_skewed If `TRUE`, factors are reorientated to ensure that the
#' loadings of each factor (i.e., the source signal matrix) have positive skew.
#' Helps ensure that the most influential features for each factor are
#' positively associated with it.
#'
#' @param scale_components If `TRUE`, the loadings are standardised (to have a
#' mean of 0 and standard deviation of 1).
#'
#' @param scale_reduced If `TRUE`, the reduced data (mixture matrix) are
#' standardised (to have a mean of 0 and standard deviation of 1).
#'
#' @param n_runs The number of times to run ICA to estimate factors and quantify
#' stability. Ignored if `use_stability` is `FALSE`.
#'
#' @param BPPARAM A class containing parameters for parallel evaluation. Uses
#' \link[BiocParallel]{SerialParam} by default, running only a single
#' ICA computation at a time. Ignored if `use_stability`
#' is `FALSE`.
#'
#' @param ... Additional arguments to be passed to
#' \link[ica]{ica}.
#'
#' @inheritParams estimateFactors
#'
#' @seealso [ica::ica()], [ReducedExperiment::estimateStability()]
#'
#' @details
#' Function performs ICA for a data matrix. If `use_stability` is `TRUE`, then
#' ICA is performed multiple times with either: i) random initialisation
#' (default); or ii) bootstrap resampling of the data (if `resample` is `TRUE`).
#'
#' Note that the seed must be set if reproducibility is needed. Specifically,
#' one can use `set.seed` prior to running standard ICA
#' (`use_stability = FALSE`) or set the `RNGseed` argument of `BPPARAM` when
#' running stabilised ICA (`use_stability = TRUE`).
#'
#' The stability-based ICA algorithm is similar to the the ICASSO approach
#' (\url{https://www.cs.helsinki.fi/u/ahyvarin/papers/Himberg03.pd}) that is
#' implemented in the stabilized-ica Python package
#' (\url{https://github.com/ncaptier/stabilized-ica/tree/master}).
#'
#' In short, the stability-based algorithm consists of:
#' \itemize{
#'  \item Running ICA multiple times with either random initialisation or
#'  bootstrap resampling of the input data.
#'  \item Clustering the resulting factors across all runs based on the
#'  signature matrix.
#'  \item Calculating intra- (aics) and extra- (aecs) cluster
#'  stability, and defining the final cluster stability as `aics - aecs`.
#'  \item Calculating the cluster centrotype as the factor with the highest
#'  intra-cluster stability.
#'  \item Optionally removing factors below a specified stability threshold
#' (`stability_threshold`).
#' }
#'
#' Results from this function should be broadly similar to those generated by
#' other implementations of stabilised ICA, although they will not be identical.
#' Notable differences include:
#' \describe{
#'  \item{ICA algorithm}{Differences in the underlying implementation of
#'  ICA.}
#'  \item{Stability threshold}{The `stability_threshold` argument, if
#'  specified, removes unstable components. Such a threshold is not
#'  used by stabilized-ica.}
#'  \item{Mixture matrix recovery}{ICA is generally formulated as
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
#' @returns A list containing the following:
#' \describe{
#'  \item{M}{The mixture matrix (reduced data) with samples as rows and columns
#'  as factors.}
#'  \item{S}{The source signal matrix (loadings) with rows as features and
#'  columns as factors.}
#'  \item{stab}{If `use_stability` is TRUE, "stab" will be a component of the
#'  list. It is a vector indicating the relative stability, as described
#'  above.}
#' }
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 100 rows (features)
#' # and 20 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(100, 20, "feature", "obs")
#'
#' # Run standard ICA on the data with 5 components
#' set.seed(1)
#' ica_res <- runICA(X, nc = 5, use_stability = FALSE)
#'
#' # Run stabilised ICA on the data with 5 components (low runs for example)
#' ica_res_stab <- runICA(X, nc = 5, use_stability = TRUE, n_runs = 5,
#'                         BPPARAM = BiocParallel::SerialParam(RNGseed = 1))
#'
#' @import ica
#' @export
runICA <- function(
    X, nc,
    use_stability = FALSE, resample = FALSE,
    method = "fast",
    stability_threshold = NULL,
    center_X = TRUE, scale_X = FALSE,
    reorient_skewed = TRUE,
    scale_components = TRUE, scale_reduced = TRUE,
    n_runs = 30,
    BPPARAM = BiocParallel::SerialParam(),
    ...
) {
    if (center_X | scale_X) {
        X <- t(scale(t(X), center = center_X, scale = scale_X))
    }

    if (use_stability) {
        ica_res <- .stabilityICA(X,
            nc = nc, resample = resample,
            method = method, stability_threshold = stability_threshold,
            n_runs = n_runs, BPPARAM = BPPARAM, ...
        )
    } else {
        if (resample) {
            stop("Cannot use resampling approach when `use_stability` is FALSE")
        }

        if (!is.null(stability_threshold)) {
            stop(
                "Cannot apply `stability_threshold` when `use_stability` ",
                "is FALSE"
            )
        }

        ica_res <- list(S = ica::ica(X,
            nc = nc, method = method,
            center = FALSE, ...
        )$S)
    }

    # Reorient and scale factors before recalculating M
    if (reorient_skewed) ica_res$S <- .reorientFactors(ica_res$S)
    if (scale_components) ica_res$S <- scale(ica_res$S)
    ica_res$M <- .projectICA(X, ica_res$S)
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
#' Function for running stabilised ICA. See \link[ReducedExperiment]{runICA}.
#'
#' @import ica
#' @import BiocParallel
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.stabilityICA <- function(
    X, nc, resample, method, n_runs, BPPARAM,
    stability_threshold, BPOPTIONS = bpoptions(), return_centrotypes = TRUE,
    ...
) {
    # Run stabilized ICA in parallel (depending on BPPARAM)
    S_all <- BiocParallel::bplapply(seq_len(n_runs), .icaRandom,
        BPPARAM = BPPARAM, BPOPTIONS = BPOPTIONS, X_mat = X, nc = nc,
        method = method, resample = resample, ...
    )
    S_all <- do.call(cbind, S_all)

    # Get correlations between factors and resulting clusters
    S_cor <- abs(stats::cor(S_all))
    S_clust <- stats::cutree(stats::hclust(stats::as.dist(1 - S_cor)), k = nc)
    S_clust <- factor(S_clust)
    names(S_clust) <- colnames(S_all)

    stabilities <- c()
    centrotypes <- data.frame(matrix(
        nrow = nrow(S_all), ncol = nc,
        dimnames = list(rownames(S_all), seq_len(nc))
    ))

    for (comp in seq_len(nc)) {
        cluster_labels <- names(S_clust)[S_clust == comp]
        non_cluster_labels <-
            names(S_clust)[!names(S_clust) %in% cluster_labels]

        # Average intra-cluster similarity - average extra-cluster similarity
        aics <- mean(S_cor[cluster_labels, cluster_labels])
        aecs <- mean(S_cor[cluster_labels, non_cluster_labels])
        stabilities <- c(stabilities, aics - aecs)

        which_is_centrotype <-
            which.max(apply(S_cor[cluster_labels, cluster_labels], 2, sum))
        centrotypes[[comp]] <- S_all[, cluster_labels[which_is_centrotype]]
    }

    if (!return_centrotypes)
        return(list(
            stab = stabilities, S_all = S_all, S_clust = S_clust, S_cor = S_cor
        ))

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

#' Parallelisable function for running a stabilized ICA iteration
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.icaRandom <- function(X_mat, i, nc, method, resample, ...) {
    # Randomly initialises ICA
    Rmat <- matrix(stats::rnorm(nc**2), nrow = nc, ncol = nc)

    if (resample) {
        X_mat <- X_mat[, sample(ncol(X_mat), replace = TRUE)]
    }

    # Get ICA loadings for given initialisation
    # (and possibly bootstrap resample)
    S <- ica::ica(
        X_mat,
        nc = nc,
        method = method,
        center = FALSE,
        Rmat = Rmat,
        ...
    )$S

    colnames(S) <- paste0("iteration_", i, "_", seq_len(ncol(S)))
    if (ncol(S) != nc) {
        warning(
            "ICA did not return expected number of factors, potentially ",
            "indicating a rank deficiency in the input"
        )
    }

    return(S)
}

#' Ensure factors have positive skew
#'
#' Reorientates factors based on their skew. Generally ensures that the most
#' aligned genes are positively associated with the factor.
#'
#' @import moments
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.reorientFactors <- function(S) {

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
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.projectICA <- function(newdata, S) {
    # M <- t(newdata) %*% t(MASS::ginv(S))
    M <- t(newdata) %*% S

    colnames(M) <- colnames(S)

    return(M)
}

#' Estimate stability of factors as a function of the number of components
#'
#' Estimates the stability of factors over a range of component numbers to
#' aid in the identification of the optimal factor number. Based on the
#' Most Stable Transcriptome Dimension (MSTD) approach (see `details`).
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
#' @param mean_stability_threshold A threshold for the mean stability of factors.
#'
#' @param verbose If `TRUE`, shows a progress bar that updates for each
#' number of components tested. Note that the time taken may not be linear,
#' because the time taken to run ICA generally increases with the number
#' of components.
#'
#' @param ... Additional arguments to be passed to
#' \link[ReducedExperiment]{runICA}.
#'
#' @inheritParams runICA
#' @inheritParams estimateFactors
#'
#' @seealso [ReducedExperiment::runICA()], [ReducedExperiment::plotStability()]
#'
#' @details
#' Runs the stability-based ICA algorithm
#' (see \link[ReducedExperiment]{runICA}) for a range of component numbers.
#' Estimates stability for each, allowing for selection of the optimal
#' number of components to be used for ICA. The results of this function
#' can be plotted by \link[ReducedExperiment]{plotStability}.
#'
#' This algorithm is based on the Most Stable Transcriptome
#' Dimension (MSTD) approach
#' (\url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4112-9}).
#'
#' The function automatically selects a number of components based on
#' `mean_stability_threshold`. However, this choice should be made after
#' visualisating the stabilities as a function of the number of components,
#' which may be done using \link[ReducedExperiment]{plotStability}. The
#' aformentioned MSTD paper provides additional context and advice for choosing
#' the number of components based on these results.
#'
#' @returns Returns a list containing:
#' \describe{
#'  \item{stability}{A data.frame indicating factor stabilities as a function
#'  of the number of components.}
#'  \item{selected_nc}{a naive estimate for the optimal number of components
#'  based on the `mean_stability_threshold`.}
#' }
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 200 rows (features)
#' # and 100 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(200, 100, "feature", "obs")
#'
#' # Estimate stability across 10 to 30 components
#' # Note: We could have provided a SummarizedExperiment object instead of a matrix
#' stab_res_1 <- estimateStability(
#'     X,
#'     min_components = 10,
#'     max_components = 30,
#'     n_runs = 5,
#'     verbose = FALSE
#' )
#'
#' @export
estimateStability <- function(
    X, min_components = 10, max_components = 60, by = 2,
    n_runs = 30, resample = FALSE, mean_stability_threshold = NULL,
    center_X = TRUE, scale_X = FALSE, assay_name = "normal",
    BPPARAM = BiocParallel::SerialParam(), verbose = TRUE,
    ...
) {
    if (inherits(X, "SummarizedExperiment")) {
        X <- assay(X, assay_name)
    }

    if (dim(X)[2] < max_components) {
        stop("Number of samples must be greater than max_components")
    }

    stabilities <- data.frame()

    if (verbose) {
        tpb <- utils::txtProgressBar(
            min = min_components, max = max_components,
            initial = min_components, style = 3
        )
    }

    for (nc in seq(from = min_components, to = max_components, by = by)) {
        ica_res <- runICA(X,
            nc = nc, center_X = center_X, scale_X = scale_X,
            use_stability = TRUE, resample = resample, BPPARAM = BPPARAM,
            method = "fast", n_runs = n_runs, ...
        )

        stabilities <- rbind(stabilities, data.frame(
            nc = nc, component_name = names(ica_res$stab),
            component_number = as.numeric(gsub(
                "factor_", "",
                names(ica_res$stab)
            )), stability = ica_res$stab
        ))

        if (verbose) utils::setTxtProgressBar(tpb, nc)
    }

    if (verbose) close(tpb)

    select_nc <- .selectNC(stabilities, mean_stability_threshold)
    return(list("stability" = stabilities, "selected_nc" = select_nc))
}

#' Estimates an appropriate number of components
#'
#' @author Jack Gisby
#'
#' @noRd
#' @keywords internal
.selectNC <- function(stabilities, mean_stability_threshold) {
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
}

#' Plot component stability as a function of the number of components
#'
#' Plots the results of \link[ReducedExperiment]{estimateStability}. See this
#' function's documentation for more information.
#'
#' @param stability The results of \link[ReducedExperiment]{estimateStability}.
#'
#' @param plot_path The path at which the plot will be saved
#'
#' @param stability_threshold Plots a stability threshold, below which
#' components can be pruned by \link[ReducedExperiment]{runICA}.
#'
#' @param mean_stability_threshold Plots a stability threshold, which is used
#' by \link[ReducedExperiment]{estimateStability} to provide a naive estimate
#' for the optimal number of components.
#'
#' @param height The height of the plot, to be passed to \link[ggplot2]{ggsave}.
#'
#' @param width The width of the plot, to be passed to \link[ggplot2]{ggsave}.
#'
#' @param ... Additional arguments to be passed to \link[ggplot2]{ggsave}.
#'
#' @returns Returns a list of three plots as `ggplot2` objects:
#' \describe{
#'  \item{combined_plot}{The two other plots combined with patchwork.}
#'  \item{stability_plot}{A plot in which each line indicates stability
#'   as a function of the number of components. A line is shown for each
#'   number of components tested.}
#'  \item{mean_plot}{The average component stability as a function of the
#'  number of components.}
#' }
#'
#' @author Jack Gisby
#'
#' @examples
#' # Get a random matrix with rnorm, with 200 rows (features)
#' # and 100 columns (observations)
#' X <- ReducedExperiment:::.makeRandomData(200, 100, "feature", "obs")
#'
#' # Estimate stability across 10 to 30 components
#' stab_res <- estimateStability(
#'     X,
#'     min_components = 10,
#'     max_components = 30,
#'     n_runs = 5,
#'     verbose = FALSE
#' )
#'
#' # Intracluster stability similar to extracluster since this is random data
#' plotStability(stab_res)$combined_plot
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
plotStability <- function(
    stability, plot_path = NULL,
    stability_threshold = NULL, mean_stability_threshold = NULL,
    height = 4, width = 10,
    ...
) {
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

    mean_stab_plot <- ggplot(
        stabilities_agg,
        aes(!!sym("nc"), !!sym("stability"), group = 1)
    ) +
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
