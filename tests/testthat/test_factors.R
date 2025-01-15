context("reduce_data")

test_that("Estimate factors", {
    # Use real data from airway package
    airway <- .getAirwayData()
    assayNames(airway) <- c("counts", "normalised")  # Test with alternative assay name
    airway_fe <- estimateFactors(airway, nc = 2, assay_name = "normalised")

    expect_equal(dim(airway_fe), c("Features" = 18870, "Samples" = 8, "Components" = 2))
    expect_equal(reduced(airway_fe)[1, ], c("factor_1" = 0.0562, "factor_2" = -1.2809), tolerance = 1e-3)
    expect_equal(loadings(airway_fe)[1, ], c("factor_1" = -0.741, "factor_2" = -0.275), tolerance = 1e-3)
    expect_equivalent(length(airway_fe@center), nFeatures(airway_fe))
})

test_that("runICA matches estimateFactors", {
    # Use real data from airway package
    airway <- .getAirwayData()

    set.seed(1)
    airway_fe <- estimateFactors(airway, nc = 2)

    set.seed(1)
    airway_ica <- runICA(assay(airway, "normal"), nc = 2)

    expect_equal(reduced(airway_fe), airway_ica$M, tolerance = 1e-3)
    expect_equal(loadings(airway_fe), airway_ica$S, tolerance = 1e-3)
})

test_that("Stability", {
    # Use real data from airway package

    set.seed(2)
    airway <- .getAirwayData(n_features = 500)
    assayNames(airway) <- c("counts", "normalised")  # Test with alternative assay name

    stability_res <- estimateStability(airway, n_runs = 5, min_components = 2, max_components = 4, by = 1, mean_stability_threshold = 0.9, assay_name = "normalised")
    stability_res_bootstrap <- estimateStability(airway, n_runs = 5, min_components = 2, max_components = 4, by = 1, mean_stability_threshold = 0.9, resample = TRUE, assay_name = "normalised")

    expect_equal(stability_res$selected_nc, 4)
    expect_equal(stability_res_bootstrap$selected_nc, NULL)

    expect_equal(colnames(stability_res$stability), c("nc", "component_name", "component_number", "stability"))
    expect_equal(colnames(stability_res_bootstrap$stability), colnames(stability_res$stability))

    plotStability(stability_res)

    airway_fe <- estimateFactors(airway, nc = 2, use_stability = TRUE, n_runs = 5, BPPARAM = BiocParallel::SerialParam(RNGseed = 1), assay_name = "normalised")

    expect_equal(reduced(airway_fe)[1, ], c("factor_1" = -1.11680, "factor_2" = 0.00221), tolerance = 1e-3)
    expect_equal(loadings(airway_fe)[1, ], c("factor_1" = -0.0196, "factor_2" = -0.6465), tolerance = 1e-3)
    expect_equal(names(stability(airway_fe)), componentNames(airway_fe))
    expect_true(all(stability(airway_fe) > 0.99))

    airway_fe_bootstrap <- estimateFactors(airway, nc = 2, use_stability = TRUE, n_runs = 5, resample = TRUE, BPPARAM = BiocParallel::SerialParam(RNGseed = 1), assay_name = "normalised")

    expect_equal(reduced(airway_fe_bootstrap)[1, ], c("factor_1" = -0.6462, "factor_2" = -1.2427), tolerance = 1e-2)
    expect_equal(loadings(airway_fe_bootstrap)[1, ], c("factor_1" = -0.6539, "factor_2" = -0.0837), tolerance = 1e-2)
    expect_equal(names(stability(airway_fe_bootstrap)), componentNames(airway_fe_bootstrap))
    expect_true(all(stability(airway_fe_bootstrap) > 0.2))
})
