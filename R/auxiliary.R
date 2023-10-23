#' Turn rnorm into a matrix
.makeRandomData <- function(r, c, rname, cname, seed=1) {
    set.seed(seed)

    m <- matrix(rnorm(n = r * c), nrow = r, ncol = c)

    rownames(m) <- as.character(paste0(rname, "_", 1:r))
    colnames(m) <- as.character(paste0(cname, "_", 1:c))

    return(m)
}
