context("Mortality")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj$q_matrix
  Q <- as.matrix(Q)

  # Extract R
  R <- diag(nrow(Q))
  diag(R) <- samc_obj@data@t_abs

  R_list <- lapply(split(samc_obj@data@c_abs, col(samc_obj@data@c_abs)),
                   function(x){
                     mat <- diag(length(x))
                     diag(mat) <- x
                     return(mat)
                   })

  R_list <- c(list(total = R), R_list)

  # Create an indentity matrix
  I <- diag(nrow(Q))

  # Fundamental matrix
  F_mat <- solve(I - Q)

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
    # Make sure absorption components add up
    samc_obj$override <- TRUE
    result <- mortality(samc_obj)
    samc_obj$override <- FALSE
    expect_equal(as.vector(result[[1]]), as.vector(Reduce('+', result) - result[[1]]))

    # Make sure all absorption components match
    base_result <- lapply(R_list, function(x) F_mat %*% x)
    mapply(function(x, y) expect_equal(as.vector(x), as.vector(y)),
           base_result, result)
  })


  test_that("Testing mortality(samc, origin)", {
    # Make sure absorption components add up
    result <- mortality(samc_obj, origin = row_vec[1])
    expect_equal(as.vector(result[[1]]), as.vector(Reduce('+', result) - result[[1]]))

    # Make sure named and unnamed results match up
    result_char <- mortality(samc_obj, origin = as.character(row_vec[1]))
    expect_equal(result[[1]], result_char[[1]])

    # Make sure all absorption components match
    base_result <- lapply(R_list, function(x) F_mat %*% x)
    mapply(function(x, y) expect_equal(as.vector(x[row_vec[1], ]), as.vector(y)),
           base_result, result)
  })

  test_that("Testing mortality(samc, dest)", {
    # Make sure absorption components add up
    result <- mortality(samc_obj, dest = col_vec[1])
    expect_equal(as.vector(result[[1]]), as.vector(Reduce('+', result) - result[[1]]))

    # Make sure named and unnamed results match up
    result_char <- mortality(samc_obj, dest = as.character(col_vec[1]))
    expect_equal(result[[1]], result_char[[1]])

    # Make sure all absorption components match
    base_result <- lapply(R_list, function(x) F_mat %*% x)
    mapply(function(x, y) expect_equal(as.vector(x[, col_vec[1]]), as.vector(y)),
           base_result, result)
  })

  test_that("Testing mortality(samc, origin, dest)", {
    # Make sure absorption components add up
    vector_result <- mortality(samc_obj, origin = row_vec, des = col_vec)
    expect_equal(as.vector(vector_result[[1]]), as.vector(Reduce('+', vector_result) - vector_result[[1]]))

    # Make sure named and unnamed results match up
    vector_result_char <- mortality(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))
    expect_equal(vector_result[[1]], vector_result_char[[1]])


    base_result <- lapply(R_list, function(x) F_mat %*% x)
    for (i in 1:length(row_vec)) {
      # Test single pair version against paired vector version
      r <- mortality(samc_obj, origin = row_vec[i], dest = col_vec[i])
      mapply(function(x, y) expect_equal(x[i], y),
             vector_result, r)

      # Test against base result
      mapply(function(x, y) expect_equal(x, y[row_vec[i], col_vec[i]], check.names = FALSE),
             r, base_result)
    }
  })

  test_that("Testing mortality(samc, occ)", {
    # Make sure absorption components add up
    result <- mortality(samc_obj, occ = test$occ)
    expect_equal(as.vector(result[[1]]), as.vector(Reduce('+', result) - result[[1]]))

    # Make sure all absorption components match
    base_result <- lapply(R_list, function(x) pv %*% F_mat %*% x)
    mapply(function(x, y) expect_equal(as.vector(x), as.vector(y)),
           base_result, result)
    })
}
