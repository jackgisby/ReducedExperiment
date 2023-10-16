context("ReducedSet")

dir.create("tempTestOutput")

test_that("Build ReducedSet", {

  r <- matrix(nrow = 2)
  r[1] <- 5
  rownames(r) <- c("1", "2")
  colnames(r) <- "1"

  l <- matrix()
  rownames(l) <- "1"
  colnames(l) <- "1"

  e <- matrix(ncol = 2)
  rownames(e) <- "1"
  colnames(e) <- c("1", "2")

  x <- new("ReducedSet", exprs=e, reducedData=r, L=l)
  x

  expect_equal(dim(x), c("Features" = 1, "Samples" = 2))
})

unlink("tempTestOutput", recursive = TRUE)
