for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj$q_matrix
  Q <- as.matrix(Q)

  # Prepare the occupancy data
  pv <- as_pv(test$init)

  # Calculate P^t
  Pt <- mat_power(Q, time)


  # Run the tests
  test_that("Testing distribution(samc, time)", {
    samc_obj$override <- TRUE
    result <- distribution(samc_obj, time = time)
    samc_obj$override <- FALSE

    base_result <- Pt

    expect_equal(dim(result), dim(base_result))
    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that("Testing distribution(samc, origin, time)", {
    result <- distribution(samc_obj, origin = row_vec[1], time = time)
    result_char <- distribution(samc_obj, origin = as.character(row_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- Pt[row_vec[1], ]

    expect_equal(result, base_result, check.names = FALSE)
  })

  test_that("Testing distribution(samc, origin, time_vec)", {
    result <- distribution(samc_obj, origin = row_vec[1], time = time_vec)
    result_char <- distribution(samc_obj, origin = as.character(row_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- mat_power(Q, time_vec[i])[row_vec[1], ]

      expect_equal(result[[i]], base_result, check.names = FALSE)
    }
  })

  test_that("Testing distribution(samc, dest, time)", {
    result <- distribution(samc_obj, dest = col_vec[1], time = time)
    result_char <- distribution(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- Pt[, col_vec[1]]

    expect_equal(result, base_result, check.names = FALSE)
  })

  test_that("Testing distribution(samc, dest, time_vec)", {
    result <- distribution(samc_obj, dest = col_vec[1], time = time_vec)
    result_char <- distribution(samc_obj, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- mat_power(Q, time_vec[i])[, col_vec[1]]

      expect_equal(result[[i]], base_result, check.names = FALSE)
    }
  })

  test_that("Testing distribution(samc, origin, dest, time)", {
    result <- distribution(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time)
    result_char <- distribution(samc_obj, origin = as.character(row_vec[1]), dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- Pt[row_vec[1], col_vec[1]]

    expect_equal(result, base_result)
  })

  test_that("Testing distribution(samc, origin, dest, time_vec)", {
    result <- distribution(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time_vec)
    result_char <- distribution(samc_obj, origin = as.character(row_vec[1]), dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- mat_power(Q, time_vec[i])[row_vec[1], col_vec[1]]

      expect_equal(result[[i]], base_result)
    }
  })

  test_that("Testing distribution(samc, init, time)", {
    result <- distribution(samc_obj, init = test$init, time = time)

    base_result <- pv %*% (Pt)

    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that("Testing distribution(samc, init, time_vec)", {
    result <- distribution(samc_obj, init = test$init, time = time_vec)

    for (i in 1:length(time_vec)) {
      base_result <- pv %*% mat_power(Q, time_vec[i])

      expect_equal(result[[i]], as.vector(base_result))
    }
  })
}
