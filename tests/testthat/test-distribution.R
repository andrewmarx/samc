context("Visitation")


for(test in testlist) {
  # Create the samc object
  samc_obj <- test$samc

  # Extract Q
  Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
  Q <- as.matrix(Q)

  # Prepare the occupancy data
  occ_ras <- raster::raster(test$occ)
  pv <- as.vector(occ_ras)
  pv <- pv[is.finite(pv)]

  # Calculate P^t
  Pt <- Q
  for (i in 2:time) {
    Pt <- Pt %*% Q
  }


  # Run the tests
  test_that("Testing distribution(samc, time)", {

    r1 <- distribution(samc_obj, time = time)

    r2 <- Pt

    # Verify
    expect_equal(dim(r1), dim(r2))
    expect_equal(as.vector(r1), as.vector(r2))
  })

  test_that("Testing distribution(samc, origin, time)", {

    r1 <- distribution(samc_obj, origin = row_vec[1], time = time)

    r2 <- Pt[row_vec[1], ]

    # Verify
    expect_equal(r1, r2)
  })

  test_that("Testing distribution(samc, origin, time_vec)", {

    r1 <- distribution(samc_obj, origin = row_vec[1], time = time_vec)

    for (i in 1:length(time_vec)) {
      pt <- Q
      for (j in 2:time_vec[i]) {
        pt <- pt %*% Q
      }
      r2 <- pt[row_vec[1], ]

      # Verify
      expect_equal(r1[[i]], r2)
    }
  })

  test_that("Testing distribution(samc, dest, time)", {

    r1 <- distribution(samc_obj, dest = col_vec[1], time = time)

    r2 <- Pt[, col_vec[1]]

    # Verify
    expect_equal(r1, r2)
  })

  test_that("Testing distribution(samc, dest, time_vec)", {

    r1 <- distribution(samc_obj, dest = col_vec[1], time = time_vec)

    for (i in 1:length(time_vec)) {
      pt <- Q
      for (j in 2:time_vec[i]) {
        pt <- Q %*% pt
      }
      r2 <- pt[, col_vec[1]]

      # Verify
      expect_equal(r1[[i]], r2)
    }
  })

  test_that("Testing distribution(samc, origin, dest, time)", {

    r1 <- distribution(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time)

    r2 <- Pt[row_vec[1], col_vec[1]]

    # Verify
    expect_equal(r1, r2)
  })

  test_that("Testing distribution(samc, origin, dest, time_vec)", {

    r1 <- distribution(samc_obj, origin = row_vec[1], dest = col_vec[1], time = time_vec)

    for (i in 1:length(time_vec)) {
      pt <- Q
      for (j in 2:time_vec[i]) {
        pt <- pt %*% Q
      }
      r2 <- pt[row_vec[1], col_vec[1]]

      # Verify
      expect_equal(r1[[i]], r2)
    }
  })

  test_that("Testing distribution(samc, occ, time)", {
    r1 <- distribution(samc_obj, occ = test$occ, time = time)

    r2 <- pv %*% (Pt)

    # Verify
    expect_equal(as.vector(r1), as.vector(r2))
  })

  test_that("Testing distribution(samc, occ, time_vec)", {
    r1 <- distribution(samc_obj, occ = test$occ, time = time_vec)

    for (i in 1:length(time_vec)) {
      pt <- Q
      for (j in 2:time_vec[i]) {
        pt <- pt %*% Q
      }
      r2 <- pv %*% pt

      # Verify
      expect_equal(r1[[i]], as.vector(r2))
    }
  })
}
