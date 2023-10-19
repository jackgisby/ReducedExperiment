context("ReducedSet")

dir.create("tempTestOutput")

#' i (features), j (samples), k components)
.createRandomisedReducedSet <- function(i, j, k) {
    return(new(
        "ReducedSet",
        exprs = .makeRandomData(i, j, "gene", "sample"),
        reduced = .makeRandomData(j, k, "sample", "factor")
    ))
}

test_that("Build ReducedSet", {

    i <- 300
    j <- 100
    k <- 10

    rrs <- .createRandomisedReducedSet(i=i, j=j, k=k)

    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(exprs(rrs)), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))

    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features"=6, "Samples"=41, "Components"=2))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(paste0("sample_", 50:90), sampleNames(rrs_subset))
})

unlink("tempTestOutput", recursive = TRUE)
