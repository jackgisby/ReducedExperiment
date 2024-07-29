context("reduce_data")

test_that("Estimate factors", {
    # Use real data from airway package
    airway <- .get_airway_data()

    airway_fe <- estimate_factors(airway, nc = 2, seed = 1)

    expect_equal(dim(airway_fe), c("Features" = 18870, "Samples" = 8, "Components" = 2))
    expect_equal(reduced(airway_fe)[1, ], c("factor_1" = 0.0562, "factor_2" = -1.2809), tolerance = 1e-3)
    expect_equal(loadings(airway_fe)[1, ], c("factor_1" = -0.741, "factor_2" = -0.275), tolerance = 1e-3)
    expect_equivalent(length(airway_fe@center), nFeatures(airway_fe))
})

test_that("run_ica matches estimate_factors", {
    # Use real data from airway package
    airway <- .get_airway_data()

    airway_fe <- estimate_factors(airway, nc = 2, seed = 1)
    airway_ica <- run_ica(assay(airway, "normal"), nc = 2, seed = 1)

    expect_equal(reduced(airway_fe), airway_ica$M, tolerance = 1e-3)
    expect_equal(loadings(airway_fe), airway_ica$S, tolerance = 1e-3)
})

test_that("Stability ICA", {
    # Use real data from airway package
    airway <- .get_airway_data(n_features = 500)

    airway_fe <- estimate_factors(airway, nc = 2, seed = 1, use_stability = TRUE, n_runs = 5)

    expect_equal(reduced(airway_fe)[1, ], c("factor_1" = -1.11680, "factor_2" = 0.00221), tolerance = 1e-3)
    expect_equal(loadings(airway_fe)[1, ], c("factor_1" = -0.0196, "factor_2" = -0.6465), tolerance = 1e-3)
    expect_equal(names(stability(airway_fe)), componentNames(airway_fe))
    expect_true(all(stability(airway_fe) > 0.99))

    airway_fe_bootstrap <- estimate_factors(airway, nc = 2, seed = 1, use_stability = TRUE, n_runs = 5, resample = TRUE)

    expect_equal(reduced(airway_fe_bootstrap)[1, ], c("factor_1" = -1.38, "factor_2" = -1.56), tolerance = 1e-2)
    expect_equal(loadings(airway_fe_bootstrap)[1, ], c("factor_1" = -0.2433, "factor_2" = 0.0809), tolerance = 1e-2)
    expect_equal(names(stability(airway_fe_bootstrap)), componentNames(airway_fe_bootstrap))
    expect_true(all(stability(airway_fe_bootstrap) > 0.2))
})

test_that("Estimate stability", {
    # Use real data from airway package
    airway <- .get_airway_data(n_features = 500)

    stability_res <- estimate_stability(airway, seed = 1, n_runs = 5, min_components = 2, max_components = 4, by = 1, mean_stability_threshold = 0.9)
    stability_res_bootstrap <- estimate_stability(airway, seed = 1, n_runs = 5, min_components = 2, max_components = 4, by = 1, mean_stability_threshold = 0.9, resample = TRUE)

    expect_equal(stability_res$selected_nc, 4)
    expect_equal(stability_res_bootstrap$selected_nc, NULL)

    expect_equal(colnames(stability_res$stability), c("nc", "component_name", "component_number", "stability"))
    expect_equal(colnames(stability_res_bootstrap$stability), colnames(stability_res$stability))

    plot_stability(stability_res)
})
