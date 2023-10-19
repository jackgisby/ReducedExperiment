context("reduce_data")

dir.create("tempTestOutput")

test_that("Reduce data", {
    factor_set <- reduce_data(.makeRandomData(50, 30, "gene", "sample"), method="ICA", nc=5)

    expect_equal(dim(factor_set), c("Features" = 50, "Samples" = 30, "Components" = 5))
})

unlink("tempTestOutput", recursive = TRUE)
