context("Conditional Passage Time")

br_function <- function(samc, col) {
  # Calculate the results based on De Sanctis and de Koning 2018
  Q <- samc$q_matrix
  Q <- as.matrix(Q)

  qj <- Q[-col, col]
  Qj <- Q[-col, -col]

  I <- diag(nrow(Qj))

  r <- rowSums(samc$r_matrix)
  r <- r[-col]

  R <- cbind(r, qj)

  f <- solve(I - Qj)

  b <- as.matrix(f %*% R)
  bdg <- Matrix::sparseMatrix(i = 1:nrow(b),
                              j = 1:nrow(b),
                              x = b[, 2],
                              index1 = TRUE)

  bdg <- as.matrix(bdg)

  res <- solve(bdg) %*% f %*% bdg %*% rep(1, nrow(bdg))
  return(as.numeric(res))
}

for(test in testlist) {
  # TODO cond_passage does not work in all cases yet. Remove this when it does
  if (!(test$id %in% c(1, 2))) next

  # Create the samc object
  samc_obj <- test$samc

  # Create a version from P matrix
  samc_p <- samc(samc_obj$p_matrix)

  # Run the tests
  test_that("Testing cond_passage(samc, dest)", {

    base_result <- br_function(samc_obj, col_vec[1])

    r1 <- cond_passage(samc_p, dest = col_vec[1])
    r2 <- cond_passage(samc_p, dest = as.character(col_vec[1]))

    # Verify
    expect_equal(dim(r1), dim(base_result))
    expect_equal(as.vector(r1), as.vector(base_result))
    expect_equal(r1, r2)
  })

  test_that("Testing cond_passage(samc, origin, dest)", {

    vector_result <- cond_passage(samc_p, row_vec, col_vec)
    vector_result_char <- cond_passage(samc_p, as.character(row_vec), as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      base_result <- br_function(samc_obj, col_vec[i])

      r <- cond_passage(samc_p, origin = row_vec[i], dest = col_vec[i])

      # Offset needed because a column is removed and sometimes the row # is greater than or equal to that col #
      if (row_vec[i] > col_vec[i]) {
        offset = -1
      } else if (row_vec[i] == col_vec[i]) {
        offset = length(base_result) # forces the vector lookup out of bounds to produce a NA
      } else {
        offset = 0
      }

      expect_equal(r, base_result[row_vec[!!i] + offset])
      expect_equal(vector_result[i], r)
    }
  })
}
