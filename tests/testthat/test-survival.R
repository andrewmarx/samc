for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Fundamental matrix
  F_mat <- ref_f(samc_obj)

  # Prepare the occupancy data
  pv <- as_pv(test$init)

  ones <- rep(1, nrow(F_mat))


  # Run the tests
  test_that(sprintf("Testing survival(samc) [scenario %d]", test$id), {

    r1 <- survival(samc_obj)

    r2 <- F_mat %*% ones

    # Verify equality
    expect_equal(as.vector(r1), as.vector(r2))
  })

  test_that(sprintf("Testing survival(samc, init) [scenario %d]", test$id), {
    # Calculate psi*z using survival(samc, init)
    r1 <- survival(samc_obj, test$init)

    r2 <- pv %*% F_mat %*% ones

    # Verify equality
    expect_equal(as.vector(r1), as.vector(r2))
  })
}
