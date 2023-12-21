context("ReducedExperiment")

test_that("Build and subset", {

    i <- 300
    j <- 100
    k <- 10

    rrs <- .createRandomisedReducedExperiment(i=i, j=j, k=k)

    # Test reduced experiment
    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))

    # Subset and re-test
    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features" = 6, "Samples" = 41, "Components" = 2))
    expect_equal(nComponents(rrs_subset), c("Components" = 2))
    expect_equal(nSamples(rrs_subset), c("Samples" = 41))
    expect_equal(nFeatures(rrs_subset), c("Features" = 6))
    expect_true(ncol(rrs_subset) == nSamples(rrs_subset))
    expect_true(nrow(rrs_subset) == nFeatures(rrs_subset))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(paste0("sample_", 50:90), sampleNames(rrs_subset))

    # Now test an empty object
    rrs_empy <- ReducedExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0))
})

test_that("Access and replace reduced data", {

    i <- 300
    j <- 100
    k <- 10

    reduced_data <- .makeRandomData(j, k, "sample", "factor")

    # Make a ReducedExperiment and save original reduced data to an object
    rrs <- ReducedExperiment(
        assays = list("normal"=.makeRandomData(i, j, "gene", "sample")),
        reduced = reduced_data
    )

    # Expect that the reduced method returns the original data
    expect_equal(reduced(rrs), reduced_data)

    # Replacing a value should work
    reduced(rrs)[3, 5] <- 5
    expect_equal(reduced(rrs)[3, 5], 5)

    reduced(rrs)[3, 5] <- reduced_data[3, 5]
    expect_equal(reduced(rrs), reduced_data)

    # This should work
    reduced(rrs) <- reduced_data[, 1:5]

    # This should not work (validity should fail because different number of samples)
    expect_error((function() {reduced(rrs) <- reduced_data[1:5 ,]})())
})

test_that("Access and replace component names", {

    rrs <- .createRandomisedReducedExperiment(i=300, j=100, k=10)

    expect_equal(componentNames(rrs), paste0("factor_", 1:10))
    expect_equal(colnames(reduced(rrs)), paste0("factor_", 1:10))

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
})

test_that("Access and replace sample names", {

    rrs <- .createRandomisedReducedExperiment(i=300, j=100, k=10)

    expect_equal(sampleNames(rrs), paste0("sample_", 1:100))
    expect_equal(rownames(colData(rrs)), paste0("sample_", 1:100))
    expect_equal(colnames(assay(rrs, 1)), paste0("sample_", 1:100))
    expect_equal(rownames(reduced(rrs)), paste0("sample_", 1:100))

    sampleNames(rrs)[5] <- "new_name"
    expect_equal(sampleNames(rrs)[5], "new_name")
    expect_equal(rownames(colData(rrs))[5], "new_name")
    expect_equal(colnames(assay(rrs, 1))[5], "new_name")
    expect_equal(rownames(reduced(rrs))[5], "new_name")
})

test_that("Access and replace feature names", {
    rrs <- .createRandomisedReducedExperiment(i=300, j=100, k=10)

    expect_equal(featureNames(rrs), paste0("gene_", 1:300))
    expect_equal(rownames(assay(rrs, 1)), paste0("gene_", 1:300))

    featureNames(rrs)[5] <- "new_name"
    expect_equal(featureNames(rrs)[5], "new_name")
    expect_equal(rownames(assay(rrs, 1))[5], "new_name")
})
