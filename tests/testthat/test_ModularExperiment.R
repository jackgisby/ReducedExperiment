context("ModularExperiment")

test_that("Build and subset", {

    i <- 300
    j <- 100
    k <- 10

    rrs <- .createRandomisedModularExperiment(i=i, j=j, k=k)

    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))

    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features" = 6, "Samples" = 41, "Components" = 2))
    expect_equal(nComponents(rrs_subset), c("Components" = 2))
    expect_equal(nComponents(rrs_subset), nModules(rrs_subset))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))

    expect_equal(paste0("sample_", 50:90), sampleNames(rrs_subset))

    rrs_empy <- ModularExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0))
    expect_equal(assignments(rrs_empy), character())
})

test_that("Access and replace assignments", {

})

test_that("Access and replace component/module names", {

})

test_that("Eigengene calculation / projection / prediction", {

})
