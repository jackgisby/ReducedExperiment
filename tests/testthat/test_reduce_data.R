context("reduce_data")

dir.create("tempTestOutput")

test_that("Reduce data", {
    scale_X <- TRUE

    random_expr <- .makeRandomData(50, 30, "gene", "sample")
    factor_exp <- run_ica(random_expr, method="ICA", nc=5, scale_X=scale_X)

    expect_equal(dim(factor_exp), c("Features" = 50, "Samples" = 30, "Components" = 5))

    reprojected_data <- project_data(factor_exp, random_expr)
    expect_equal(reduced(factor_exp), reprojected_data)

    reprojected_data <- project_data(factor_exp, t(scale(t(random_expr), scale=scale_X)), new_data_is_transformed=TRUE)
    expect_equal(reduced(factor_exp), reprojected_data)

    # project_data(factor_exp, .makeRandomData(50, 20, "gene", "sample"))
})

unlink("tempTestOutput", recursive = TRUE)
