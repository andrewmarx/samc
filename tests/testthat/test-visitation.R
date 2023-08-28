context("Visitation")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj$q_matrix
  Q <- as.matrix(Q)

  # Prepare the occupancy data
  occ_ras = raster::raster(test$init)
  pv = as.vector(occ_ras)
  pv = pv[is.finite(pv)]


  #
  # Short-term tests
  #

  Qt <- diag(nrow(Q))
  base_result = Qt
  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    base_result <- base_result + Qt
  }

  test_that("Testing visitation(samc, time)", {
    samc_obj$override = TRUE
    r = visitation(samc_obj, time = time)
    samc_obj$override = FALSE

    expect_equal(dim(r), dim(base_result))
    expect_equal(as.vector(r), as.vector(base_result))
  })

  test_that("Testing visitation(samc, origin, time)", {
    result = visitation(samc_obj, origin = row_vec[1], time = time)
    result_char = visitation(samc_obj, origin = as.character(row_vec[1]), time = time)
    expect_equal(result, result_char)

    expect_equal(as.vector(result), as.vector(base_result[row_vec[1], ]))
  })

  test_that("Testing visitation(samc, dest, time)", {
    result = visitation(samc_obj, dest = col_vec[1], time = time)
    result_char = visitation(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    expect_equal(as.vector(result), as.vector(base_result[, col_vec[1]]))
  })

  test_that("Testing visitation(samc, init, time)", {
    result = visitation(samc_obj, init = test$init, time = time)

    r = pv %*% base_result
    expect_equal(as.vector(result), as.vector(r))
  })


  #
  # Long-term tests
  #

  I = diag(nrow(Q))
  base_result = solve(I - Q)

  test_that("Testing visitation(samc)", {
    samc_obj$override <- TRUE
    r <- visitation(samc_obj)
    samc_obj$override <- FALSE

    expect_equal(dim(r), dim(base_result))
    expect_equal(as.vector(r), as.vector(base_result))
  })

  test_that("Testing visitation(samc, origin)", {
    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i])
      r_char <- visitation(samc_obj, origin = as.character(row_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, base_result[row_vec[i], ], check.names = FALSE)
    }
  })

  test_that("Testing visitation(samc, dest)", {
    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, dest = col_vec[i])
      r_char <- visitation(samc_obj, dest = as.character(col_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, base_result[, col_vec[i]], check.names = FALSE)
    }
  })

  test_that("Testing visitation(samc, origin, dest)", {
    # TODO fix multi-location inputs
    skip("Skip dues to multi-location input")
    vector_result <- visitation(samc_obj, origin = row_vec, dest = col_vec)
    vector_result_char <- visitation(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in 1:length(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, base_result[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })


  # TODO visitation_net() tests
}
