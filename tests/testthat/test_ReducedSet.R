context("ReducedSet")

dir.create("tempTestOutput")

test_that("TRUE is TRUE", {
  expect_equal(TRUE, TRUE)
})

unlink("tempTestOutput", recursive = TRUE)
