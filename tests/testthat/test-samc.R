context("SAMC test")


test_that("samc p matrix row sums equal 1 without fidelity data", {

  # Create the samc object and get the row sums of the p matrix
  samc_obj <- samc(testlist[[1]]$res, testlist[[1]]$abs, tr_fun = function(x) 1/mean(x))
  rs <- Matrix::rowSums(samc_obj@p)

  # Create a vector of the expected result based on the number of non-NA cells
  # in the original data and an extra entry for the last row of the p matrix
  length <- sum(!is.na(as.vector(testlist[[1]]$res)))
  v <- numeric(length + 1)
  v[] <- 1

  # Verify equality
  expect_equal(rs, v)
})

test_that("samc p matrix row sums equal 1 with fidelity data", {

  # Create the samc object and get the row sums of the p matrix
  samc_obj <- samc(testlist[[1]]$res, testlist[[1]]$abs, testlist[[1]]$fid, tr_fun = function(x) 1/mean(x))
  rs <- Matrix::rowSums(samc_obj@p)

  # Create a vector of the expected result based on the number of non-NA cells
  # in the original data and an extra entry for the last row of the p matrix
  length <- sum(!is.na(as.vector(testlist[[1]]$res)))
  v <- numeric(length + 1)
  v[] <- 1

  # Verify equality
  expect_equal(rs, v)
})
