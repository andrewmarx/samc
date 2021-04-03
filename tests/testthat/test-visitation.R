context("Visitation")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
  Q <- as.matrix(Q)

  # Create an indentity matrix
  I <- diag(nrow(Q))

  base_result <- solve(I - Q)

  # Run the tests
  test_that("Testing visitation(samc)", {
    r <- visitation(samc_obj)

    expect_equal(dim(r), dim(base_result))
    expect_equal(as.vector(r), as.vector(base_result))
  })

  test_that("Testing visitation(samc, origin)", {
    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i])
      r_char <- visitation(samc_obj, origin = as.character(row_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, base_result[row_vec[i], ], check.names = FALSE)
    }
  })

  test_that("Testing visitation(samc, dest)", {
    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, dest = col_vec[i])
      r_char <- visitation(samc_obj, dest = as.character(col_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, base_result[, col_vec[i]], check.names = FALSE)
    }
  })

  test_that("Testing visitation(samc, origin, dest)", {
    vector_result <- visitation(samc_obj, origin = row_vec, dest = col_vec)
    vector_result_char <- visitation(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })
}
