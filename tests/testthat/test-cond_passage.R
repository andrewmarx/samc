context("Conditional Passage Time")


for(test in testlist) {
  # TODO cond_passage does not work in all cases yet. Remove this when it does
  if (!(test$id %in% c(1, 2))) next

  # Create the samc object
  samc_obj <- test$samc

  # Create a version from P matrix
  samc_p <- samc(p_mat = samc_obj@p)


  # Calculate the results based on De Sanctis and de Koning 2018
  Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
  Q <- as.matrix(Q)

  qj <- Q[-col_vec[1], col_vec[1]]
  Qj <- Q[-col_vec[1], -col_vec[1]]

  I <- diag(nrow(Qj))

  r <- samc_obj@p[-nrow(samc_obj@p), ncol(samc_obj@p)]
  r <- r[-col_vec[1]]

  R <- cbind(r, qj)

  f <- solve(I - Qj)

  b <- as.matrix(f %*% R)
  bdg <- Matrix::sparseMatrix(i = 1:nrow(b),
                              j = 1:nrow(b),
                              x = b[, 2],
                              index1 = TRUE)

  bdg <- as.matrix(bdg)

  result <- solve(bdg) %*% f %*% bdg %*% rep(1, nrow(bdg))
  result <- as.numeric(result)

  # Run the tests
  test_that("Testing cond_passage(samc, dest)", {

    r1 <- cond_passage(samc_p, dest = col_vec[1])

    # Verify
    expect_equal(dim(r1), dim(result))
    expect_equal(as.vector(r1), as.vector(result))
  })

  test_that("Testing cond_passage(samc, origin, dest)", {

    r1 <- cond_passage(samc_p, origin = row_vec[1], dest = col_vec[1])

    r_vec <- cond_passage(samc_p, row_vec, col_vec)

    # Verify
    expect_equal(as.vector(r1), as.vector(result[row_vec[1]]))
    for (i in 1:length(row_vec)) {
      r1 <- cond_passage(samc_p, row_vec[i], col_vec[i])
      expect_equal(r_vec[i], r1)
    }
  })
}
