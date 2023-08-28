context("Conditional Passage Time")

br_function <- function(samc, col) {
  # Calculate the results based on De Sanctis and de Koning 2018
  Q <- samc$q_matrix
  Q <- as.matrix(Q)

  qj <- Q[-col, col]
  Qj <- Q[-col, -col]

  I <- diag(nrow(Qj))

  r <- samc@data@t_abs
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

test_num = 0
for(test in testlist) {
  test_num = test_num + 1

  # TODO cond_passage does not work in all cases yet. Remove this when it does
  if (!(test$id %in% c(1, 2))) next

  # Create the samc object
  samc_obj <- test$samc

  # Create a version from P matrix
  samc_p <- samc(samc_obj$p_matrix)
  samc_p@clumps = samc_obj@clumps # TODO: remove when creation of samc object from P matrix calculates clumps


  # Run the tests
  test_that(paste("Testing cond_passage(samc, dest):", test_num), {

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

  test_that(paste("Testing cond_passage(samc, origin, dest)", test_num), {
    # TODO fix multi-location inputs
    skip("Skip dues to multi-location input")
    vector_result <- cond_passage(samc_p, origin = row_vec, dest = col_vec)
    vector_result_char <- cond_passage(samc_p, origin = as.character(row_vec), dest = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      base_result <- cond_passage(samc_obj, dest = col_vec[i])

      r <- cond_passage(samc_p, origin = row_vec[i], dest = col_vec[i])

      expect_equal(r, unname(base_result[row_vec[i]]))
      expect_equal(vector_result[i], r)
    }
  })
}
