for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- ref_q(samc_obj)

  # Prepare the occupancy data
  pv <- as_pv(test$init)


  #
  # Short-term tests
  #

  st_base <- sum_powers(Q, time)

  test_that(sprintf("Testing visitation(samc, time) [scenario %d]", test$id), {
    samc_obj$override <- TRUE
    r <- visitation(samc_obj, time = time)
    samc_obj$override <- FALSE

    expect_equal(dim(r), dim(st_base))
    expect_equal(as.vector(r), as.vector(st_base))
  })

  test_that(sprintf("Testing visitation(samc, origin, time) [scenario %d]", test$id), {
    result <- visitation(samc_obj, origin = row_vec[1], time = time)
    result_char <- visitation(samc_obj, origin = as.character(row_vec[1]), time = time)
    expect_equal(result, result_char)

    expect_equal(as.vector(result), as.vector(st_base[row_vec[1], ]))
  })

  test_that(sprintf("Testing visitation(samc, dest, time) [scenario %d]", test$id), {
    result <- visitation(samc_obj, dest = col_vec[1], time = time)
    result_char <- visitation(samc_obj, dest = as.character(col_vec[1]), time = time)
    expect_equal(result, result_char)

    expect_equal(as.vector(result), as.vector(st_base[, col_vec[1]]))
  })

  test_that(sprintf("Testing visitation(samc, init, time) [scenario %d]", test$id), {
    result <- visitation(samc_obj, init = test$init, time = time)

    r <- pv %*% st_base
    expect_equal(as.vector(result), as.vector(r))
  })


  #
  # Long-term tests
  #

  lt_base <- ref_f(samc_obj)

  test_that(sprintf("Testing visitation(samc) [scenario %d]", test$id), {
    samc_obj$override <- TRUE
    r <- visitation(samc_obj)
    samc_obj$override <- FALSE

    expect_equal(dim(r), dim(lt_base))
    expect_equal(as.vector(r), as.vector(lt_base))
  })

  test_that(sprintf("Testing visitation(samc, origin) [scenario %d]", test$id), {
    for (i in seq_along(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i])
      r_char <- visitation(samc_obj, origin = as.character(row_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, lt_base[row_vec[i], ], check.names = FALSE)
    }
  })

  test_that(sprintf("Testing visitation(samc, dest) [scenario %d]", test$id), {
    for (i in seq_along(row_vec)) {
      r <- visitation(samc_obj, dest = col_vec[i])
      r_char <- visitation(samc_obj, dest = as.character(col_vec[i]))

      expect_equal(r, r_char)
      expect_equal(r, lt_base[, col_vec[i]], check.names = FALSE)
    }
  })

  test_that(sprintf("Testing visitation(samc, origin, dest) [scenario %d]", test$id), {
    vector_result <- visitation(samc_obj, origin = row_vec, dest = col_vec)
    vector_result_char <- visitation(samc_obj, origin = as.character(row_vec), dest = as.character(col_vec))

    expect_equal(vector_result, vector_result_char)

    for (i in seq_along(row_vec)) {
      r <- visitation(samc_obj, origin = row_vec[i], dest = col_vec[i])

      expect_equal(vector_result[i], r)
      expect_equal(r, lt_base[row_vec[i], col_vec[i]], check.names = FALSE)
    }
  })


  # TODO visitation_net() tests
}
