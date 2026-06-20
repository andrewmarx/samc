for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- ref_q(samc_obj)

  # Create an indentity matrix
  I <- diag(nrow(Q))

  # Prepare the occupancy data
  pv <- as_pv(test$init)

  # Pre-calc
  f <- ref_f(samc_obj)
  fdg <- I
  diag(fdg) <- 1/diag(f)


  #Run the tests
  test_that(sprintf("Testing dispersal(samc, dest, time) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, dest = col_vec[1], time = time)
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]

    base_result <- sum_powers(Qj, time) %*% qj

    expect_equal(as.vector(result)[-col_vec[1]], as.vector(base_result))
  })

  test_that(sprintf("Testing dispersal(samc, dest, time_vec) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, dest = col_vec[1], time = time_vec)
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]
    for (i in seq_along(time_vec)) {
      base_result <- sum_powers(Qj, time_vec[i]) %*% qj

      expect_equal((result[[i]])[-col_vec[1]], as.vector(base_result))
    }
  })

  test_that(sprintf("Testing dispersal(samc, init, dest, time) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, init = test$init, dest = col_vec[1], time = time)
    result_char <- dispersal(samc_obj, init = test$init, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]

    base_result <- pv[-col_vec[1]] %*% (sum_powers(Qj, time) %*% qj)

    expect_equal(result, as.numeric(base_result))
  })

  test_that(sprintf("Testing dispersal(samc, init, dest, time_vec) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, init = test$init, dest = col_vec[1], time = time_vec)
    result_char <- dispersal(samc_obj, init = test$init, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    for (i in seq_along(time_vec)) {
      Qj <- Q[-col_vec[1],-col_vec[1]]

      base_result <- pv[-col_vec[1]] %*% (sum_powers(Qj, time_vec[i]) %*% qj)

      expect_equal(result[[i]], as.numeric(base_result))
    }
  })

  test_that(sprintf("Testing dispersal(samc) [scenario %d]", test$id), {
    samc_obj$override <- TRUE
    result <- dispersal(samc_obj)
    samc_obj$override <- FALSE

    base_result <- (f - I) %*% fdg

    expect_equal(dim(result), dim(base_result))
    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that(sprintf("Testing dispersal(samc, origin) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, origin = row_vec[1])
    samc_obj@.cache$dgf_exists <- FALSE
    samc_obj$threads <- 2
    result_par <- dispersal(samc_obj, origin = row_vec[1])
    samc_obj@.cache$dgf_exists <- FALSE
    samc_obj$threads <- 1
    expect_equal(result, result_par)
    result_char <- dispersal(samc_obj, origin = as.character(row_vec[1]))
    expect_equal(result, result_char)

    base_result <- (f - I) %*% fdg

    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], ]))
  })

  test_that(sprintf("Testing dispersal(samc, dest) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, dest = col_vec[1])
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]))
    expect_equal(result, result_char)

    base_result <- (f - I) %*% fdg

    # Verify
    expect_equal(as.vector(result), as.vector(base_result[, col_vec[1]]))
  })

  test_that(sprintf("Testing dispersal(samc, origin, dest) [scenario %d]", test$id), {
    base_result <- (f - I) %*% fdg

    vector_result <- dispersal(samc_obj, origin = row_vec, dest = col_vec)
    vector_result_char <- dispersal(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))
    expect_equal(vector_result, vector_result_char)

    for (i in seq_along(row_vec)) {
      r <- dispersal(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })

  test_that(sprintf("Testing dispersal(samc, init) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, init = test$init)
    samc_obj@.cache$dgf_exists <- FALSE
    samc_obj$threads <- 2
    result_par <- dispersal(samc_obj, init = test$init)
    samc_obj@.cache$dgf_exists <- FALSE
    samc_obj$threads <- 1

    expect_equal(result, result_par)

    base_result <- pv %*% (f - I) %*% fdg

    # Verify
    expect_equal(as.vector(result), as.vector(base_result))
  })


  test_that(sprintf("Testing dispersal(samc, init, dest) [scenario %d]", test$id), {
    result <- dispersal(samc_obj, init = test$init, dest = col_vec[1])
    result_char <- dispersal(samc_obj, init = test$init, dest = as.character(col_vec[1]))
    expect_equal(result, result_char)

    base_result <- pv %*% (f - I) %*% fdg

    # Verify
    expect_equal(result, as.vector(base_result)[col_vec[1]])
  })
}
