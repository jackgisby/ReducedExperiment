context("samples")

test_that("Associate components", {
    airway <- .get_airway_data(n_features = 500)
    airway_fe <- estimate_factors(airway, nc = 2, seed = 1)

    fe_res <- associate_components(airway_fe, "~ cell + dex")

    expect_true(all(paste0("factor_", 1:2) %in% fe_res$anovas$component))
    expect_true(all(paste0("factor_", 1:2) %in% fe_res$summaries$component))
    expect_true(all(paste0("factor_", 1:2) %in% names(fe_res$models)))
    expect_true(inherits(fe_res$models[[1]], "lm"))

    airway_me <- identify_modules(airway, verbose = 0, powers = 21)
    me_res <- associate_components(airway_me, "~ dex + (1|cell)", method = "lmer")

    expect_true(all(paste0("module_", 0:5) %in% me_res$anovas$component))
    expect_true(all(paste0("module_", 0:5) %in% me_res$summaries$component))
    expect_true(all(paste0("module_", 0:5) %in% names(me_res$models)))
    expect_true(inherits(me_res$models[[1]], "lmerModLmerTest"))
})
