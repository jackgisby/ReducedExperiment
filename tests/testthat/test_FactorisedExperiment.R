context("FactorisedExperiment")

test_that("Build and subset", {

    i <- 300
    j <- 100
    k <- 10

    # Create random FactorisedExperiment
    rrs <- .createRandomisedFactorisedExperiment(i=i, j=j, k=k)

    # Check dimensions and names are as expected
    expect_equal(dim(rrs), c("Features" = i, "Samples" = j, "Components" = k))
    expect_equal(dim(rrs), c(nFeatures(rrs), nSamples(rrs), nComponents(rrs)))

    expect_equal(colnames(assay(rrs, "normal")), sampleNames(rrs))
    expect_equal(rownames(reduced(rrs)), sampleNames(rrs))

    # Do the same after slicing
    rrs_subset <- rrs[5:10, 50:90, 1:2]
    expect_equal(dim(rrs_subset), c("Features" = 6, "Samples" = 41, "Components" = 2))
    expect_equal(nComponents(rrs_subset), c("Components" = 2))
    expect_equal(sampleNames(rrs_subset), paste0("sample_", 50:90))
    expect_equal(rownames(reduced(rrs_subset)), sampleNames(rrs_subset))
    expect_equal(featureNames(rrs_subset), paste0("gene_", 5:10))
    expect_equal(rownames(loadings(rrs_subset)), featureNames(rrs_subset))

    # Do the same with an empty experiment
    rrs_empy <- FactorisedExperiment()
    expect_equal(dim(rrs_empy), c("Features" = 0, "Samples" = 0, "Components" = 0))
    expect_equal(loadings(rrs_empy), matrix(0, 0, 0))
    expect_equal(reduced(rrs_empy), matrix(0, 0, 0), check.attributes = FALSE)
    expect_equal(stability(rrs_empy), NULL)

})

test_that("Access and replace component names", {

    rrs <- .createRandomisedFactorisedExperiment(i=300, j=100, k=10)

    expect_equal(componentNames(rrs), paste0("factor_", 1:10))
    expect_equal(colnames(reduced(rrs)), paste0("factor_", 1:10))
    expect_equal(colnames(loadings(rrs)), paste0("factor_", 1:10))

    componentNames(rrs)[5] <- "new_name"
    expect_equal(componentNames(rrs)[5], "new_name")
    expect_equal(colnames(reduced(rrs))[5], "new_name")
    expect_equal(colnames(loadings(rrs))[5], "new_name")
})

test_that("Access and replace feature names", {

    rrs <- .createRandomisedFactorisedExperiment(i=300, j=100, k=10)

    expect_equal(featureNames(rrs), paste0("gene_", 1:300))
    expect_equal(rownames(assay(rrs, 1)), paste0("gene_", 1:300))
    expect_equal(rownames(loadings(rrs)), paste0("gene_", 1:300))

    featureNames(rrs)[5] <- "new_name"
    expect_equal(featureNames(rrs)[5], "new_name")
    expect_equal(rownames(assay(rrs, 1))[5], "new_name")
    expect_equal(rownames(loadings(rrs))[5], "new_name")
})

test_that("Access and replace loadings", {

    i <- 300
    j <- 100
    k <- 10

    loadings_data <- .makeRandomData(i, k, "gene", "factor")

    # Make a FactorisedExperiment and save original loadings data to an object
    rrs <- FactorisedExperiment(
        assays = list("normal" = .makeRandomData(i, j, "gene", "sample")),
        reduced = .makeRandomData(j, k, "sample", "factor"),
        loadings = loadings_data,
        stability = NULL
    )

    # Expect that the loadings method returns the original data
    expect_equal(loadings(rrs), loadings_data)

    # Replacing a value should work
    loadings(rrs)[3, 5] <- 5
    expect_equal(loadings(rrs)[3, 5], 5)

    loadings(rrs)[3, 5] <- loadings_data[3, 5]
    expect_equal(loadings(rrs), loadings_data)

    # This should not work (validity should fail because there are a different number of factors in loadings vs. reduced slots)
    expect_error((function() {loadings(rrs) <- loadings_data[, 1:5]})())

    # Neither should this (validity should fail because different number of samples)
    expect_error((function() {loadings(rrs) <- loadings_data[1:5 ,]})())
})

test_that("Access and replace stability", {

    rrs <- .createRandomisedFactorisedExperiment(i=300, j=100, k=10)

    expect_equal(stability(rrs), NULL)

    # Can change stability given that it matches with expected factor length/names
    stability(rrs) <- 1:10
    stability(rrs) <- setNames(1:10, paste0("factor_", 1:10))

    # Should not be able to create vector with non-matching length/names
    expect_error((function() {stability(rrs) <- setNames(1:10, paste0("module_", 1:10))})())
    expect_error((function() {stability(rrs) <- setNames(1:5, paste0("factor_", 1:5))})())
    expect_error((function() {stability(rrs) <- 1:5})())
})

test_that("Predict and project", {

    # Use real data from airway package
    airway <- .get_airway_data()
    airway_fe <- estimate_factors(airway, nc=2, seed=1, scale_components=FALSE, reorient_skewed=FALSE)

    # Check that projecting the data reproduces the original results
    for (input_type in c("se", "matrix", "data.frame")) {

        if (input_type == "se") {
            newdata <- airway_fe
        } else if (input_type == "matrix") {
            newdata <- as.matrix(assay(airway, "normal"))
        } else if (input_type == "data.frame") {
            newdata <- as.data.frame(assay(airway, "normal"))
        }

        for (projection_function in c(projectData, predict)) {

            res <- projection_function(airway_fe, newdata)

            if (input_type == "se") res <- reduced(res)

            expect_equal(res, reduced(airway_fe), check.attributes = FALSE)
        }
    }

    # Check that projection method is equivalent to original results from ica::ica
    set.seed(1)
    ica_res <- ica::ica(t(scale(t(assay(airway_fe, "normal")), center=TRUE, scale=FALSE)),
                        nc=2, method="fast", center=FALSE)

    expect_equal(scale(ica_res$M), scale(matrix(reduced(airway_fe), ncol = 2)), check.attributes = FALSE)
})

test_that("Get aligned features", {

    rrs <- .createRandomisedFactorisedExperiment(i=300, j=100, k=10)

    # Check z cutoff
    aligned_features <- getAlignedFeatures(rrs, z_cutoff=3, n_features=NULL,
                                           feature_id_col="rownames", format="list")

    expect_equal(aligned_features, list(
        "factor_2"="gene_195",
        "factor_4"="gene_265",
        "factor_5"=c("gene_95", "gene_242"),
        "factor_9"="gene_60"
    ))

    # Get top three factors for each
    aligned_features <- getAlignedFeatures(rrs, z_cutoff=3, n_features=3,
                                           feature_id_col="rownames", format="list")

    expect_equal(aligned_features, list(
        "factor_1"=c("gene_232", "gene_274", "gene_206"),
        "factor_2"=c("gene_195", "gene_146", "gene_61"),
        "factor_3"=c("gene_243", "gene_56", "gene_15"),
        "factor_4"=c("gene_265", "gene_75", "gene_63"),
        "factor_5"=c("gene_95", "gene_242", "gene_248"),
        "factor_6"=c("gene_237", "gene_33", "gene_99"),
        "factor_7"=c("gene_41", "gene_52", "gene_83"),
        "factor_8"=c("gene_32", "gene_267", "gene_104"),
        "factor_9"=c("gene_60", "gene_139", "gene_206"),
        "factor_10"=c("gene_112", "gene_26", "gene_293")
    ))

    # Should work the same way for these data
    expect_equal(getAlignedFeatures(rrs, z_cutoff=3, n_features=3),
                 getAlignedFeatures(rrs, z_cutoff=NULL, n_features=3))

    # Sanity check of factor 1 z cutoff approach
    aligned_features <- getAlignedFeatures(rrs, z_cutoff=1.5, n_features=NULL)
    expect_equal(aligned_features[[1]], names(loadings(rrs, scale_loadings=TRUE)[,1][abs(loadings(rrs, scale_loadings=TRUE)[,1]) > 1.5]))
})

test_that("Get gene IDs", {

    airway <- .get_airway_data(n_features=500)
    airway_fe <- estimate_factors(airway, nc=2, seed=1, use_stability=FALSE, method="imax")

    airway_fe <- getGeneIDs(airway_fe)

    expect_true("hgnc_symbol" %in% colnames(rowData(airway_fe)))
    expect_true("entrezgene_id" %in% colnames(rowData(airway_fe)))
    expect_true(mean(is.na(rowData(airway_fe)$hgnc_symbol)) < 0.05)
    expect_true(mean(is.na(rowData(airway_fe)$entrezgene_id)) < 0.3)
})
