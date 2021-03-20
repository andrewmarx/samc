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
    expect_equal(r1, r2[row_vec[1], ], check.names = FALSE)
  })

  test_that("Testing visitation(samc, dest)", {

    r1 <- visitation(samc_obj, dest = col_vec[1])

    r2 <- solve(I - Q)

    # Verify equality
    expect_equal(r1, r2[, col_vec[1]], check.names = FALSE)
  })

  test_that("Testing visitation(samc, origin, dest)", {

    base_result <- solve(I - Q)
    vector_result <- visitation(samc_obj, origin = row_vec, des = col_vec)

    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })
}
