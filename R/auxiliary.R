#' Turn rnorm into a matrix
.makeRandomData <- function(r, c, rname, cname, seed=1) {
    set.seed(seed)

    m <- matrix(rnorm(n = r * c), nrow = r, ncol = c)

    rownames(m) <- as.character(paste0(rname, "_", 1:r))
    colnames(m) <- as.character(paste0(cname, "_", 1:c))

    return(m)
}

#' i (features), j (samples), k components)
.createRandomisedReducedExperiment <- function(i, j, k, seed=1) {
    return(ReducedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample", seed=seed)),
        reduced = .makeRandomData(j, k, "sample", "factor", seed=seed)
    ))
}

#' i (features), j (samples), k components)
.createRandomisedFactorisedExperiment <- function(i, j, k, seed=1) {
    return(FactorisedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample", seed=seed)),
        reduced = .makeRandomData(j, k, "sample", "factor", seed=seed),
        loadings = .makeRandomData(i, k, "gene", "factor", seed=seed)
    ))
}

#' i (features), j (samples), k components)
.createRandomisedModularExperiment <- function(i, j, k, seed=1) {

    assignments <- paste0("gene_", 1:i)
    names(assignments) <- paste0("module_", round(runif(i, 1, k), 0))

    return(ModularExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample", seed=seed)),
        reduced = .makeRandomData(j, k, "sample", "module", seed=seed),
        assignments = assignments,
        loadings = .makeRandomData(i, 1, "gene", "gene", seed=seed)[,1]
    ))
}

.get_airway_data <- function(n_features=NULL) {

    # Get data
    data(airway, package="airway")

    # Remove genes that aren't expressed
    airway <- airway[apply(assay(airway, "counts"), 1, function(x) {all(x != 0)}) ,]

    # Remove genes at random for faster tests
    if (!is.null(n_features)) {
        set.seed(2)
        airway <- airway[sample(nrow(airway), n_features) ,]
    }

    # Do basic log transformation
    assay(airway, "normal") <- log(assay(airway, "counts") + 0.1)

    return(airway)
}
