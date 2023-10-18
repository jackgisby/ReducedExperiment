context("reduce_data")

dir.create("tempTestOutput")

test_that("Reduce data", {
    reduce_data(.makeRandomData(i, j, "gene", "sample"), method = "ICA")
})

unlink("tempTestOutput", recursive = TRUE)
