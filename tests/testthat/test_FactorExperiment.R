context("FactorisedExperiment")

dir.create("tempTestOutput")

#' i (features), j (samples), k components)
.createRandomisedFactorisedExperiment <- function(i, j, k) {
    return(FactorisedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        loadings = .makeRandomData(i, k, "gene", "factor")
    ))
}

test_that("Build FactorisedExperiment", {

    i <- 300
    j <- 100
    k <- 10

    rrs <- .createRandomisedFactorisedExperiment(i=i, j=j, k=k)

    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))

    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features"=6, "Samples"=41, "Components"=2))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(paste0("sample_", 50:90), sampleNames(rrs_subset))

    rrs_empy <- FactorisedExperiment()
})

unlink("tempTestOutput", recursive = TRUE)
