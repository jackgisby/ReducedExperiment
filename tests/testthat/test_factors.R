context("reduce_data")

test_that("Reduce data", {

    random_expr <- .makeRandomData(50, 30, "gene", "sample")
    factor_exp <- run_ica(random_expr, method="ICA", nc=5, scale_X=FALSE)

    expect_equal(dim(factor_exp), c("Features" = 50, "Samples" = 30, "Components" = 5))

    reprojected_data <- project_data(factor_exp, random_expr)
    expect_equal(reduced(factor_exp), reprojected_data)

    reprojected_data <- project_data(factor_exp, t(scale(t(random_expr), scale=FALSE)), new_data_is_centered=TRUE)
    expect_equal(reduced(factor_exp), reprojected_data)

    # project_data(factor_exp, .makeRandomData(50, 20, "gene", "sample"))
})
