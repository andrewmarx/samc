br_function <- function(samc, col) {
  # Calculate the results based on De Sanctis and de Koning 2018
  Q <- ref_q(samc)

  qj <- Q[-col, col]
  Qj <- Q[-col, -col]

  I <- diag(nrow(Qj))

  r <- samc@data@t_abs
  r <- r[-col]

  R <- cbind(r, qj)

  f <- solve(I - Qj)

  b <- as.matrix(f %*% R)
  bdg <- Matrix::sparseMatrix(i = seq_len(nrow(b)),
                              j = seq_len(nrow(b)),
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
  samc_p@clumps = samc_obj@clumps # TODO: remove when creation of samc object from P matrix calculates clumps


  # Run the tests
  test_that(sprintf("Testing cond_passage(samc, dest) [scenario %d]", test$id), {

    base_result <- br_function(samc_obj, col_vec[1])

    r1 <- cond_passage(samc_p, dest = col_vec[1])
    r2 <- cond_passage(samc_p, dest = as.character(col_vec[1]))

    r1 <- r1[-col_vec[1]]
    r2 <- r2[-col_vec[1]]

    # Verify
    expect_equal(dim(r1), dim(base_result))
    expect_equal(as.vector(r1), as.vector(base_result))
    expect_equal(r1, r2)
  })

  test_that(sprintf("Testing cond_passage(samc, origin, dest) [scenario %d]", test$id), {
    vector_result <- cond_passage(samc_p, origin = row_vec, dest = col_vec)
    vector_result_char <- cond_passage(samc_p, origin = as.character(row_vec), dest = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in seq_along(row_vec)) {
      base_result <- cond_passage(samc_obj, dest = col_vec[i])

      r <- cond_passage(samc_p, origin = row_vec[i], dest = col_vec[i])

      expect_equal(r, unname(base_result[row_vec[i]]))
      expect_equal(vector_result[i], r)
    }
  })
}
