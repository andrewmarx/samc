context("Mortality")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj$q_matrix
  Q <- as.matrix(Q)

  # Extract R
  R <- diag(nrow(Q))
  diag(R) <- rowSums(samc_obj$r_matrix)

  # Create an indentity matrix
  I <- diag(nrow(Q))

  # Prepare the occupancy data
  occ_ras <- raster::raster(test$occ)
  pv <- as.vector(occ_ras)
  pv <- pv[is.finite(pv)]


  # Run the tests
  test_that("Testing mortality(samc, time)", {
    samc_obj$override <- TRUE
    result <- mortality(samc_obj, time = time)
    samc_obj$override <- FALSE

    base_result <- diag(nrow(Q))

    Qt <- diag(nrow(Q))

    for (i in 1:(time - 1)) {
      Qt <- Qt %*% Q
      base_result <- base_result + Qt
    }

    base_result <- base_result %*% R

    # Verify
    expect_equal(dim(result), dim(base_result))
    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that("Testing mortality(samc, origin, time)", {
    result <- mortality(samc_obj, origin = row_vec[1], time = time)
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- diag(nrow(Q))

    Qt <- diag(nrow(Q))

    for (i in 1:(time - 1)) {
      Qt <- Qt %*% Q
      base_result <- base_result + Qt
    }

    base_result <- base_result %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], ]))
  })

  test_that("Testing mortality(samc, origin, time_vec)", {
    result <- mortality(samc_obj, origin = row_vec[1], time = time_vec)
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- diag(nrow(Q))

      Qt <- diag(nrow(Q))

      for (j in 1:(time_vec[i] - 1)) {
        Qt <- Qt %*% Q
        base_result <- base_result + Qt
      }

      base_result <- base_result %*% R

      # Verify
      expect_equal(result[[i]], as.vector(base_result[row_vec[1], ]))
    }
  })

  test_that("Testing mortality(samc, dest, time)", {
    result <- mortality(samc_obj, dest = col_vec[1], time = time)
    result_char <- mortality(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- diag(nrow(Q))

    Qt <- diag(nrow(Q))

    for (i in 1:(time - 1)) {
      Qt <- Qt %*% Q
      base_result <- base_result + Qt
    }

    base_result <- base_result %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result[, col_vec[1]]))
  })

  test_that("Testing mortality(samc, dest, time_vec)", {
    result <- mortality(samc_obj, dest = col_vec[1], time = time_vec)
    result_char <- mortality(samc_obj, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- diag(nrow(Q))

      Qt <- diag(nrow(Q))

      for (j in 1:(time_vec[i] - 1)) {
        Qt <- Qt %*% Q
        base_result <- base_result + Qt
      }

      base_result <- base_result %*% R

      # Verify
      expect_equal(result[[i]], as.vector(base_result[, col_vec[1]]))
    }
  })

  test_that("Testing mortality(samc, origin, dest, time)", {
    result <- mortality(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time)
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]), dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    base_result <- diag(nrow(Q))

    Qt <- diag(nrow(Q))

    for (i in 1:(time - 1)) {
      Qt <- Qt %*% Q
      base_result <- base_result + Qt
    }

    base_result <- base_result %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], col_vec[1]]))
  })

  test_that("Testing mortality(samc, origin, dest, time_vec)", {
    result <- mortality(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time_vec)
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]), dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    for (i in 1:length(time_vec)) {
      base_result <- diag(nrow(Q))

      Qt <- diag(nrow(Q))

      for (j in 1:(time_vec[i] - 1)) {
        Qt <- Qt %*% Q
        base_result <- base_result + Qt
      }

      base_result <- base_result %*% R

      # Verify
      expect_equal(result[[i]], as.vector(base_result[row_vec[1], col_vec[1]]))
    }
  })

  test_that("Testing mortality(samc, occ, time)", {
    result <- mortality(samc_obj, occ = test$occ, time = time)

    base_result <- I

    Qt <- diag(nrow(Q))

    for (i in 1:(time - 1)) {
      Qt <- Qt %*% Q
      base_result <- base_result + Qt
    }

    base_result <- pv %*% base_result %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that("Testing mortality(samc, occ, time_vec)", {
    result <- mortality(samc_obj, occ = test$occ, time = time_vec)

    for (i in 1:length(time_vec)) {
      base_result <- I

      Qt <- diag(nrow(Q))

      for (j in 1:(time_vec[i] - 1)) {
        Qt <- Qt %*% Q
        base_result <- base_result + Qt
      }

      base_result <- pv %*% base_result %*% R

      # Verify
      expect_equal(result[[i]], as.vector(base_result))
    }
  })

  test_that("Testing mortality(samc)", {
    samc_obj$override <- TRUE
    result <- mortality(samc_obj)
    samc_obj$override <- FALSE

    base_result <- solve(I - Q) %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result))
  })


  test_that("Testing mortality(samc, origin)", {
    base_result <- solve(I - Q) %*% R

    result <- mortality(samc_obj, origin = row_vec[1])
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]))

    expect_equal(result, result_char)
    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], ]))

    # Version for mult absorption
    # expect_equal(result$total, result_char$total)
    # expect_equal(as.vector(result$total), as.vector(base_result[row_vec[1], ]))
  })

  test_that("Testing mortality(samc, dest)", {
    base_result <- solve(I - Q) %*% R

    result <- mortality(samc_obj, dest = col_vec[1])
    result_char <- mortality(samc_obj, dest = as.character(col_vec[1]))

    expect_equal(result, result_char)
    expect_equal(as.vector(result), as.vector(base_result[, col_vec[1]]))
  })

  test_that("Testing mortality(samc, origin, dest)", {
    base_result <- solve(I - Q) %*% R

    vector_result <- mortality(samc_obj, origin = row_vec, des = col_vec)
    vector_result_char <- mortality(samc_obj, origin = as.character(row_vec), des = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      r <- mortality(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })

  test_that("Testing mortality(samc, occ)", {
    result <- mortality(samc_obj, occ = test$occ)

    base_result <- pv %*% solve(I - Q) %*% R

    # Verify
    expect_equal(as.vector(result), as.vector(base_result))
  })
}
