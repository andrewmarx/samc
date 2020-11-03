context("Visitation")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
  Q <- as.matrix(Q)

  # Create an indentity matrix
  I <- diag(nrow(Q))


  # Run the tests
  test_that("Testing visitation(samc)", {
    r1 <- visitation(samc_obj)

    r2 <- solve(I - Q)

    # Verify equality
    expect_equal(dim(r1), dim(r2))
    expect_equal(as.vector(r1), as.vector(r2))
  })

  test_that("Testing visitation(samc, origin)", {

    r1 <- visitation(samc_obj, origin = row_vec[1])

    r2 <- solve(I - Q)

    # Verify equality
    expect_equal(r1, r2[row_vec[1], ])
  })

  test_that("Testing visitation(samc, dest)", {

    r1 <- visitation(samc_obj, dest = col_vec[1])

    r2 <- solve(I - Q)

    # Verify equality
    expect_equal(r1, r2[, col_vec[1]])
  })

  test_that("Testing visitation(samc, origin, dest)", {

    r1 <- visitation(samc_obj, origin = row_vec[1], dest = col_vec[1])

    r2 <- solve(I - Q)

    # Verify equality
    expect_equal(r1, r2[row_vec[1], col_vec[1]])
  })
}
