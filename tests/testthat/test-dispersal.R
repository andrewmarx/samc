context("Dispersal")

for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj$q_matrix
  Q <- as.matrix(Q)

  # Extract R
  R <- diag(nrow(Q))
  diag(R) <- samc_obj@data@t_abs

  # Create an indentity matrix
  I <- diag(nrow(Q))

  # Prepare the occupancy data
  occ_ras <- raster::raster(test$occ)
  pv <- as.vector(occ_ras)
  pv <- pv[is.finite(pv)]

  # Pre-calc
  f <- solve(I - Q)
  fdg <- I
  diag(fdg) <- 1/diag(f)


  #Run the tests
  test_that("Testing dispersal(samc, dest, time)", {
    result <- dispersal(samc_obj, dest = col_vec[1], time = time)
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]

    Qji <- diag(nrow(Qj))
    base_result <- Qji

    for (i in 1:(time - 1)) {
      Qji <- Qji %*% Qj
      base_result <- base_result + Qji
    }

    base_result <- base_result %*% qj

    expect_equal(as.vector(result)[-col_vec[1]], as.vector(base_result))
  })

  test_that("Testing dispersal(samc, dest, time_vec)", {
    result <- dispersal(samc_obj, dest = col_vec[1], time = time_vec)
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]
    for (i in 1:length(time_vec)) {
      Qji <- diag(nrow(Qj))
      base_result <- Qji

      for (j in 1:(time_vec[i] - 1)) {
        Qji <- Qji %*% Qj
        base_result <- base_result + Qji
      }

      base_result <- base_result %*% qj

      expect_equal((result[[i]])[-col_vec[1]], as.vector(base_result))
    }
  })

  test_that("Testing dispersal(samc, occ, dest, time)", {
    result <- dispersal(samc_obj, occ = test$occ, dest = col_vec[1], time = time)
    result_char <- dispersal(samc_obj, occ = test$occ, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    Qj <- Q[-col_vec[1],-col_vec[1]]

    Qji <- diag(nrow(Qj))
    base_result <- Qji

    for (i in 1:(time - 1)) {
      Qji <- Qji %*% Qj
      base_result <- base_result + Qji
    }

    base_result <- pv[-col_vec[1]] %*% (base_result %*% qj)

    expect_equal(result, as.numeric(base_result))
  })

  test_that("Testing dispersal(samc, occ, dest, time_vec)", {
    result <- dispersal(samc_obj, occ = test$occ, dest = col_vec[1], time = time_vec)
    result_char <- dispersal(samc_obj, occ = test$occ, dest = as.character(col_vec[1]), time = time_vec)
    expect_equal(result, result_char)

    qj <- Q[-col_vec[1], col_vec[1]]

    for (i in 1:length(time_vec)) {
      Qj <- Q[-col_vec[1],-col_vec[1]]

      Qji <- diag(nrow(Qj))
      base_result <- Qji

      for (j in 1:(time_vec[i] - 1)) {
        Qji <- Qji %*% Qj
        base_result <- base_result + Qji
      }

      base_result <- pv[-col_vec[1]] %*% (base_result %*% qj)

      expect_equal(result[[i]], as.numeric(base_result))
    }
  })

  test_that("Testing dispersal(samc)", {
    samc_obj$override <- TRUE
    result <- dispersal(samc_obj)
    samc_obj$override <- FALSE

    base_result <- (f - I) %*% fdg

    expect_equal(dim(result), dim(base_result))
    expect_equal(as.vector(result), as.vector(base_result))
  })

  test_that("Testing dispersal(samc, origin)", {
    result <- dispersal(samc_obj, origin = row_vec[1])
    result_char <- dispersal(samc_obj, origin = as.character(row_vec[1]))
    expect_equal(result, result_char)

    base_result <- (f - I) %*% fdg

    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], ]))
  })

  test_that("Testing dispersal(samc, dest)", {
    result <- dispersal(samc_obj, dest = col_vec[1])
    result_char <- dispersal(samc_obj, dest = as.character(col_vec[1]))
    expect_equal(result, result_char)

    base_result <- (f - I) %*% fdg

    # Verify
    expect_equal(as.vector(result), as.vector(base_result[, col_vec[1]]))
  })

  test_that("Testing dispersal(samc, origin, dest)", {
    base_result <- (f - I) %*% fdg

    vector_result <- dispersal(samc_obj, origin = row_vec, dest = col_vec)
    vector_result_char <- dispersal(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))
    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      r <- dispersal(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })

  test_that("Testing dispersal(samc, occ)", {
    result <- dispersal(samc_obj, occ = test$occ)

    base_result <- pv %*% (f - I) %*% fdg

    # Verify
    expect_equal(as.vector(result), as.vector(base_result))
  })


  test_that("Testing dispersal(samc, occ, dest)", {
    result <- dispersal(samc_obj, occ = test$occ, dest = col_vec[1])
    result_char <- dispersal(samc_obj, occ = test$occ, dest = as.character(col_vec[1]))
    expect_equal(result, result_char)

    base_result <- pv %*% (f - I) %*% fdg

    # Verify
    expect_equal(result, as.vector(base_result)[col_vec[1]])
  })
}
