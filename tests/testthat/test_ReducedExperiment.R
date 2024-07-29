context("ReducedExperiment")

test_that("Build and subset", {
    i <- 300
    j <- 100
    k <- 10

    rrs <- .createRandomisedReducedExperiment(i = i, j = j, k = k)

    # Test reduced experiment
    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))
    expect_equal(rownames(rowData(rrs)), rownames(rrs))
    show(rrs)

    rrs@scale <- setNames(1:i, featureNames(rrs))
    rrs@center <- setNames(1:i, featureNames(rrs))

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
    expect_equal(names(rrs_subset@scale), featureNames(rrs_subset))
    expect_equal(names(rrs_subset@center), featureNames(rrs_subset))

    rownames(rrs_subset) <- paste0("123_", rownames(rrs_subset))
    expect_equal(rownames(rrs_subset)[1], "123_gene_5")

    # Now test an empty object
    rrs_empy <- ReducedExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0), check.attributes = FALSE)
    expect_equal(rrs_empy@scale, TRUE)
    expect_equal(rrs_empy@center, TRUE)
})

test_that("Access and replace reduced data", {
    i <- 300
    j <- 100
    k <- 10

    reduced_data <- .makeRandomData(j, k, "sample", "factor")

    # Make a ReducedExperiment and save original reduced data to an object
    rrs <- ReducedExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = reduced_data
    )

    # Expect that the reduced method returns the original data
    expect_equal(reduced(rrs, scale_reduced = TRUE, center_reduced = TRUE), scale(reduced_data))
    expect_equal(reduced(rrs), reduced_data)

    # Replacing a value should work
    reduced(rrs)[3, 5] <- 5
    expect_equal(reduced(rrs)[3, 5], 5)

    reduced(rrs)[3, 5] <- reduced_data[3, 5]
    expect_equal(reduced(rrs), reduced_data, check.attributes = FALSE)

    # This should work
    reduced(rrs) <- reduced_data[, 1:5]

    # This should not work (validity should fail because different number of samples)
    expect_error((function() {
        reduced(rrs) <- reduced_data[1:5, ]
    })())
})

test_that("Access and replace component names", {
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    expect_equal(componentNames(rrs), paste0("factor_", 1:10))
    expect_equal(colnames(reduced(rrs)), paste0("factor_", 1:10))

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
})

test_that("Access and replace sample names", {
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

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
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    rrs@scale <- setNames(1:300, featureNames(rrs))
    rrs@center <- setNames(1:300, featureNames(rrs))

    expect_equal(featureNames(rrs), paste0("gene_", 1:300))
    expect_equal(rownames(assay(rrs, 1)), paste0("gene_", 1:300))
    expect_equal(names(rrs@center), paste0("gene_", 1:300))
    expect_equal(names(rrs@scale), paste0("gene_", 1:300))

    featureNames(rrs)[5] <- "new_name"
    expect_equal(featureNames(rrs)[5], "new_name")
    expect_equal(rownames(assay(rrs, 1))[5], "new_name")
    expect_equal(names(rrs@center)[5], "new_name")
    expect_equal(names(rrs@scale)[5], "new_name")
})

test_that("Access and replace scale/center", {
    rrs <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10)

    # This should work
    rrs@scale <- setNames(1:300, featureNames(rrs))
    rrs@center <- setNames(1:300, featureNames(rrs))

    # This should not work (mismatch with feature names/length)
    expect_error((function() {
        rrs@scale <- setNames(1:10, paste0("module_", 1:10))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@scale <- setNames(1:5, paste0("factor_", 1:5))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@scale <- 1:5
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- setNames(1:10, paste0("module_", 1:10))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- setNames(1:5, paste0("factor_", 1:5))
        validObject(rrs)
    })())
    expect_error((function() {
        rrs@center <- 1:5
        validObject(rrs)
    })())
})

test_that("Combine ReducedExperiments with cbind", {
    rrs_a <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10, seed = 1)
    rrs_b <- .createRandomisedReducedExperiment(i = 300, j = 100, k = 10, seed = 2)

    # Objects should be cbind-able due to matching names
    rrs_a_a <- cbind(rrs_a, rrs_a)
    expect_true(validObject(rrs_a_a))
    rrs_a_b <- cbind(rrs_a, rrs_b)
    expect_true(validObject(rrs_a_b))

    expect_equal(dim(rrs_a_b), c("Features" = 300, "Samples" = 200, "Components" = 10))

    # Add scaling information to rrs_b but not a
    rrs_b@scale <- 1:300
    names(rrs_b@scale) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))

    # Should fail due to non-matching scaling slots
    expect_error(cbind(rrs_a, rrs_b))

    # Add scaling information to rrs_a and it should work again
    rrs_a@scale <- rrs_b@scale
    expect_no_error(cbind(rrs_a, rrs_b))

    # Same with centering
    rrs_b@center <- 1:300
    names(rrs_b@center) <- rownames(rrs_b)
    expect_true(validObject(rrs_b))
    expect_error(cbind(rrs_a, rrs_b))

    rrs_a@center <- rrs_b@center
    expect_no_error(cbind(rrs_a, rrs_b))
})
