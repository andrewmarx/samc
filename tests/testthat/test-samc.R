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

test_that("samc objects built from P matrices do not share a cache", {

  # Two absorbing chains with the same structure but different transition
  # probabilities. Solving one must not affect the results of the other.
  make_p <- function(q) {
    p <- matrix(0, 3, 3)
    p[1, ] <- c(0, q, 1 - q)
    p[2, ] <- c(q, 0, 1 - q)
    p[3, 3] <- 1
    methods::as(methods::as(Matrix::Matrix(p, sparse = TRUE), "CsparseMatrix"), "generalMatrix")
  }

  invisible(utils::capture.output(samc_a <- samc(make_p(0.5)), samc_b <- samc(make_p(0.9))))

  expect_false(identical(samc_a@.cache, samc_b@.cache))

  # Solve a first so that a shared cache would hold its factorization
  vis_a <- visitation(samc_a, origin = 1)
  vis_b <- visitation(samc_b, origin = 1)

  expect_equal(vis_a, as.vector(solve(diag(2) - as.matrix(samc_a$q_matrix))[1, ]), tolerance = 1e-10)
  expect_equal(vis_b, as.vector(solve(diag(2) - as.matrix(samc_b$q_matrix))[1, ]), tolerance = 1e-10)
})
