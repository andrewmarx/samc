context("Visitation")

library(raster)


# Create the samc object
samc_obj <- samc(res, abs, fid, tr_fun = function(x) 1/mean(x), override = TRUE)

# Extract Q
Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
Q <- as.matrix(Q)

# Prepare the occupancy data
occ_ras <- raster(occ)
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

  r1 <- distribution(samc_obj, origin = row, time = time)

  r2 <- Pt[row, ]

  # Verify
  expect_equal(r1, r2)
})

test_that("Testing distribution(samc, dest, time)", {

  r1 <- distribution(samc_obj, dest = col, time = time)

  r2 <- Pt[, col]

  # Verify
  expect_equal(r1, r2)
})

test_that("Testing distribution(samc, origin, dest, time)", {

  r1 <- distribution(samc_obj, origin = row, dest = col, time = time)

  r2 <- Pt[row, col]

  # Verify
  expect_equal(r1, r2)
})

test_that("Testing distribution(samc, occ, time)", {
  r1 <- distribution(samc_obj, occ = occ, time = time)

  r2 <- pv %*% (Pt)

  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})
