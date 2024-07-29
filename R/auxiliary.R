#' Turn rnorm into a matrix
#'
#' Used to generate a matrix of random numbers with dimensions r (rows) * c
#' (columns)
#'
#' @param r Number of rows
#'
#' @param c Number of columns
#'
#' @param rname Name of rows (will be appended to the row number)
#'
#' @param cname Name of columns (will be appended to the column number)
#'
#' @noRd
#' @keywords internal
.makeRandomData <- function(r, c, rname, cname) {
    set.seed(seed)

    m <- matrix(stats::rnorm(n = r * c), nrow = r, ncol = c)

    rownames(m) <- as.character(paste0(rname, "_", 1:r))
    colnames(m) <- as.character(paste0(cname, "_", 1:c))

    return(m)
}

#' Create a ReducedExperiment with random data
#'
#' Creates a \link[ReducedExperiment]{ReducedExperiment}
#' with i (features), j (samples), k components)
#'
#' @param i Number of features
#'
#' @param j Number of samples
#'
#' @param k Number of dimensionally-reduced components
#'
#' @noRd
#' @keywords internal
.createRandomisedReducedExperiment <- function(i, j, k) {
    return(ReducedExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor")
    ))
}

#' Create a FactorisedExperiment with random data
#'
#' Creates a \link[ReducedExperiment]{FactorisedExperiment}
#' with i (features), j (samples), k components)
#'
#' @param i Number of features
#'
#' @param j Number of samples
#'
#' @param k Number of factors
#'
#' @noRd
#' @keywords internal
.createRandomisedFactorisedExperiment <- function(i, j, k) {
    return(FactorisedExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        loadings = .makeRandomData(i, k, "gene", "factor")
    ))
}

#' Create a ModularExperiment with random data
#'
#' Creates a \link[ReducedExperiment]{ModularExperiment}
#' with i (features), j (samples), k components)
#'
#' @param i Number of features
#'
#' @param j Number of samples
#'
#' @param k Number of modules
#'
#' @noRd
#' @keywords internal
.createRandomisedModularExperiment <- function(i, j, k) {
    assignments <- paste0("gene_", 1:i)
    names(assignments) <- paste0("module_", round(stats::runif(i, 1, k), 0))

    return(ModularExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "module"),
        assignments = assignments,
        loadings = .makeRandomData(i, 1, "gene", "gene")[, 1]
    ))
}

#' Gets data from the airway package
#'
#' Gets data from the airway package and returns it as a SummarizedExperiment
#' with counts ("counts") and log-transformed counts ("normal") in the assay
#' slots. Non-expressed genes are removed.
#'
#' @param n_features If not NULL, the number of features (genes) to be randomly
#' selected from the original airway data. Can provide a smaller dataset for
#' testing
#'
#' @noRd
#' @keywords internal
.get_airway_data <- function(n_features = NULL) {
    # Get data
    utils::data("airway", package = "airway", envir = environment())

    # Remove genes that aren't expressed
    airway <- airway[apply(assay(airway, "counts"), 1, function(x) {
        all(x != 0)
    }), ]

    # Remove genes at random for faster tests
    if (!is.null(n_features)) {
        set.seed(2)
        airway <- airway[sample(nrow(airway), n_features), ]
    }

    # Do basic log transformation
    assay(airway, "normal") <- log(assay(airway, "counts") + 0.1)

    return(airway)
}
