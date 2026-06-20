# The expected row sums: one entry per non-NA cell plus the extra absorbing row
# of the p matrix, each summing to 1.
expected_rowsums <- rep(1, sum(!is.na(as.vector(testlist[[1]]$res))) + 1)

test_that("samc p matrix row sums equal 1 without fidelity data", {

  # Create the samc object (no fidelity) and get the row sums of the p matrix
  samc_obj <- samc(testlist[[1]]$res, testlist[[1]]$abs, model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))
  rs <- Matrix::rowSums(samc_obj$p_matrix)

  # Verify equality
  expect_equal(rs, expected_rowsums, check.names = FALSE)
})

test_that("samc p matrix row sums equal 1 with fidelity data", {

  # The shared fixture is built with fidelity data, so reuse it here
  rs <- Matrix::rowSums(testlist[[1]]$samc$p_matrix)

  # Verify equality
  expect_equal(rs, expected_rowsums, check.names = FALSE)
})
