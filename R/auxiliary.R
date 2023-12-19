#' Turn rnorm into a matrix
.makeRandomData <- function(r, c, rname, cname, seed=1) {
    set.seed(seed)

    m <- matrix(rnorm(n = r * c), nrow = r, ncol = c)

    rownames(m) <- as.character(paste0(rname, "_", 1:r))
    colnames(m) <- as.character(paste0(cname, "_", 1:c))

    return(m)
}

#' i (features), j (samples), k components)
.createRandomisedReducedExperiment <- function(i, j, k) {
    return(ReducedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor")
    ))
}

#' i (features), j (samples), k components)
.createRandomisedFactorisedExperiment <- function(i, j, k) {
    return(FactorisedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        loadings = .makeRandomData(i, k, "gene", "factor")
    ))
}

#' i (features), j (samples), k components)
.createRandomisedModularExperiment <- function(i, j, k) {

    assignments <- paste0("gene_", 1:i)
    names(assignments) <- paste0("module_", round(runif(i, 1, 4), 0))

    return(ModularExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        assignments = assignments
    ))
}
