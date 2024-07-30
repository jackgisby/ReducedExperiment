context("modules")

test_that("Identify modules", {
    set.seed(2)
    airway <- .get_airway_data(n_features = 500)

    WGCNA::disableWGCNAThreads()
    airway_me <- identify_modules(airway, verbose = 0, powers = 21)

    expect_equal(dim(airway_me), c("Features" = 500, "Samples" = 8, "Components" = 6))
    expect_equal(componentNames(airway_me), paste0("module_", 0:5))
    expect_equal(reduced(airway_me)[1, 1], 0.832, tolerance = 1e-3)
    expect_true(all(paste0("module_", 0:5) %in% names(assignments(airway_me))))
    expect_equivalent(assignments(airway_me), rownames(airway_me))
})
