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
    expect_equal(sampleNames(rrs_subset), paste0("sample_", 50:90))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(featureNames(rrs_subset), paste0("gene_", 5:10))
    expect_true(all(assignments(rrs_subset) == featureNames(rrs_subset)))

    rrs_empy <- ModularExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0))
    expect_equal(assignments(rrs_empy), character())
})

test_that("Access and replace assignments", {

    rrs <- .createRandomisedModularExperiment(i=300, j=100, k=10)

    expect_true(all(assignments(rrs) == paste0("gene_", 1:300)))

    # Can change assignments given that it matches with expected features
    assignments(rrs) <- setNames(paste0("gene_", 1:300), as.character(1:300))

    # Should not be able to create vector with non-matching length/names
    expect_error((function() {assignments(rrs) <- setNames(paste0("notgene_", 1:300), as.character(1:300))})())
    expect_error((function() {assignments(rrs) <- setNames(paste0("gene_", 1:5), as.character(1:5))})())
    expect_error((function() {assignments(rrs) <- 1:300})())
})

test_that("Access and replace component/module names", {

    rrs <- .createRandomisedModularExperiment(i=300, j=100, k=10)

    expect_equal(componentNames(rrs), paste0("module_", 1:10))
    expect_equal(colnames(reduced(rrs)), paste0("module_", 1:10))

    is_module_5 <- names(assignments(rrs)) == "module_5"

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
    expect_true(all(names(assignments(rrs))[is_module_5] == "new_name"))
    expect_true(!any(names(assignments(rrs))[!is_module_5] == "new_name"))
})

test_that("Eigengene calculation / projection / prediction", {

    # Use real data from airway package
    airway <- .get_airway_data(n_features=500)
    airway_me <- identify_modules(airway, verbose=0, powers=21)

    # Check that projecting the data reproduces the original results
    for (input_type in c("se", "matrix", "data.frame")) {

        if (input_type == "se") {
            newdata <- airway_me
        } else if (input_type == "matrix") {
            newdata <- as.matrix(assay(airway, "normal"))
        } else if (input_type == "data.frame") {
            newdata <- as.data.frame(assay(airway, "normal"))
        }

        for (projection_function in c(calcEigengenes, predict)) {

            res <- projection_function(airway_me, newdata)

            if (input_type == "se") res <- reduced(res)

            expect_equal(res, reduced(airway_me))
        }
    }
})
