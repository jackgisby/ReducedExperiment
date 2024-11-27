context("ModularExperiment")

test_that("Build and subset", {
    i <- 300
    j <- 100
    k <- 10

    set.seed(1)
    rrs <- ReducedExperiment:::.createRandomisedModularExperiment(i = i, j = j, k = k)

    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))
    expect_equal(rownames(rowData(rrs)), rownames(rrs))
    show(rrs)

    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features" = 6, "Samples" = 41, "Components" = 2))
    expect_equal(nComponents(rrs_subset), c("Components" = 2))
    expect_equal(nComponents(rrs_subset), nModules(rrs_subset))
    expect_equal(sampleNames(rrs_subset), paste0("sample_", 50:90))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(featureNames(rrs_subset), rownames(rrs_subset))
    expect_equal(featureNames(rrs_subset), paste0("gene_", 5:10))
    expect_true(all(assignments(rrs_subset) == featureNames(rrs_subset)))

    rownames(rrs_subset) <- paste0("123_", rownames(rrs_subset))
    expect_equal(rownames(rrs_subset)[1], "123_gene_5")

    rrs_empy <- ModularExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0), check.attributes = FALSE)
    expect_equal(assignments(rrs_empy), character())

    # Subset with characters
    expect_equal(
        dimnames(rrs[5:10, 50:90, 1:2]),
        dimnames(rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("module_", 1:2)])
    )

    # Test subset replacement
    rrs[5:10, 50:90, 1:2] <- rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("module_", 1:2)]
    expect_true(validObject(rrs))

    rrs[paste0("gene_", 5:10), paste0("sample_", 50:90), paste0("module_", 1:2)] <- rrs[5:10, 50:90, 1:2]
    expect_true(validObject(rrs))

    rrs_subset <- rrs[5:10, 50:90, 1:2]
    rownames(rrs_subset)[3] <- "renamed_feature"
    colnames(rrs_subset)[2] <- "renamed_sample"
    componentNames(rrs_subset)[1] <- "renamed_comp"
    rownames(rrs)[7] <- "renamed_feature"
    colnames(rrs)[51] <- "renamed_sample"
    componentNames(rrs)[1] <- "renamed_comp"

    reduced(rrs_subset)["renamed_sample", "renamed_comp"] <- 10
    names(assignments(rrs_subset))[which(assignments(rrs_subset) == "renamed_feature")] <- "renamed_module"
    loadings(rrs_subset)["renamed_feature"] <- 50

    rrs[5:10, 50:90, 1:2] <- rrs_subset
    expect_true(validObject(rrs))
    expect_true(validObject(rrs_subset))
    expect_equal(reduced(rrs)["renamed_sample", "renamed_comp"], 10)
    expect_true(names(assignments(rrs))[which(assignments(rrs) == "renamed_feature")] == "renamed_module")
    expect_true(loadings(rrs)["renamed_feature"] == 50)
})

test_that("Access and replace assignments", {
    set.seed(1)
    rrs <- .createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    expect_true(all(assignments(rrs) == paste0("gene_", 1:300)))

    # Can change assignments given that it matches with expected features
    assignments(rrs) <- setNames(paste0("gene_", 1:300), as.character(1:300))

    # Should not be able to create vector with non-matching length/names
    expect_error((function() {
        assignments(rrs) <- setNames(paste0("notgene_", 1:300), as.character(1:300))
    })())
    expect_error((function() {
        assignments(rrs) <- setNames(paste0("gene_", 1:5), as.character(1:5))
    })())
    expect_error((function() {
        assignments(rrs) <- 1:300
    })())
})

test_that("Access and replace component/module names", {
    set.seed(1)
    rrs <- .createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    expect_equal(componentNames(rrs), paste0("module_", 1:10))
    expect_equal(componentNames(rrs), moduleNames(rrs))
    expect_equal(componentNames(rrs), colnames(reduced(rrs)))

    is_module_5 <- names(assignments(rrs)) == "module_5"

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
    expect_true(all(names(assignments(rrs))[is_module_5] == "new_name"))
    expect_true(!any(names(assignments(rrs))[!is_module_5] == "new_name"))

    expect_equal(moduleNames(rrs)[5], "new_name")
    moduleNames(rrs)[3] <- "new_name_2"
    expect_equal(moduleNames(rrs)[3], "new_name_2")
})

test_that("Access and replace feature names", {
    set.seed(1)
    rrs <- .createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    expect_equal(featureNames(rrs), paste0("gene_", 1:300))
    expect_equal(rownames(assay(rrs, 1)), paste0("gene_", 1:300))
    expect_equal(names(loadings(rrs)), paste0("gene_", 1:300))

    featureNames(rrs)[5] <- "new_name"
    expect_equal(featureNames(rrs)[5], "new_name")
    expect_equal(rownames(assay(rrs, 1))[5], "new_name")
    expect_equal(names(loadings(rrs))[5], "new_name")
})

test_that("Access and replace loadings", {
    i <- 300
    j <- 100
    k <- 10

    set.seed(1)
    loadings_data <- .makeRandomData(1, i, "gene", "gene")[1, ]

    assignments_data <- paste0("gene_", 1:i)
    names(assignments_data) <- paste0("module_", round(runif(i, 1, k), 0))

    # Make a ModularExperiment and save original loadings data to an object
    set.seed(1)
    rrs <- ModularExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        assignments = assignments_data,
        loadings = loadings_data
    )

    # Expect that the loadings method returns the original data
    expect_equal(loadings(rrs), loadings_data)

    # Replacing a value should work
    loadings(rrs)[3] <- 5
    expect_true(loadings(rrs)[3] == 5)

    loadings(rrs)[3] <- loadings_data[3]
    expect_equal(loadings(rrs), loadings_data)

    # This should not work (validity should fail because there are a different number of factors in loadings vs. reduced slots)
    expect_error((function() {
        loadings(rrs) <- loadings_data[, 1:5]
    })())

    # Neither should this (validity should fail because different number of samples)
    expect_error((function() {
        loadings(rrs) <- loadings_data[1:5, ]
    })())
})

test_that("Eigengene calculation / projection / prediction", {
    # Use real data from airway package

    set.seed(2)
    airway <- .get_airway_data(n_features = 500)

    WGCNA::disableWGCNAThreads()
    airway_me <- identify_modules(airway, verbose = 0, power = 21)

    # Recalculate eigengenes using WGCNA::moduleEigengenes
    eig <- WGCNA::moduleEigengenes(t(assay(airway_me, "transformed")), setNames(names(assignments(airway_me)), assignments(airway_me)))
    colnames(eig$eigengenes) <- gsub("ME", "", colnames(eig$eigengenes))

    # Ensure method equality
    expect_equal(reduced(airway_me), scale(eig$eigengenes), check.attributes = FALSE)

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
            for (project in c(FALSE, TRUE)) {
                res <- projection_function(airway_me, newdata, project = project)

                if (input_type == "se") res <- reduced(res)

                expect_equal(scale(res), scale(reduced(airway_me)), check.attributes = FALSE)
            }
        }
    }

    # Plot and access dendrogram
    expect_true(dim(airway_me)[1] == length(dendrogram(airway_me)$order))

    plotDendro(airway_me)

    dendrogram(airway_me) <- NULL
    expect_null(dendrogram(airway_me))
})

test_that("Combine ModularExperiments with cbind and rbind", {
    set.seed(1)
    rrs_a <- ReducedExperiment:::.createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    set.seed(2)
    rrs_b <- ReducedExperiment:::.createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    # Objects should be bind-able due to matching names
    rrs_a_a <- cbind(rrs_a, rrs_a)
    expect_true(validObject(rrs_a_a))
    rrs_a_a <- rbind(rrs_a, rrs_a)
    expect_true(validObject(rrs_a_a))

    # This should fail because of non-matching loadings/assignments and reduced data
    expect_error(cbind(rrs_a, rrs_b))
    expect_error(rbind(rrs_a, rrs_b))

    # Should succeed when loadings are equivalent
    loadings(rrs_b) <- loadings(rrs_a)
    assignments(rrs_b) <- assignments(rrs_a)
    expect_true(validObject(cbind(rrs_a, rrs_b)))
    expect_error(rbind(rrs_a, rrs_b))

    # Should succeed when reduced data are equivalent
    reduced(rrs_b) <- reduced(rrs_a)
    expect_true(validObject(rbind(rrs_a, rrs_b)))
})

test_that("Get module hub genes", {
    set.seed(1)
    rrs <- .createRandomisedModularExperiment(i = 300, j = 100, k = 10)

    centrality <- getCentrality(rrs)

    for (m in componentNames(rrs)) {
        expect_equal(setdiff(centrality$feature[centrality$module == m], assignments(rrs)[names(assignments(rrs)) == m]), character())
    }

    for (i in 1:nrow(centrality)) {
        expect_equal(
            centrality$r[i],
            cor(assay(rrs)[centrality$feature[i], ], data.frame(reduced(rrs))[[centrality$module[i]]])
        )
    }
})
