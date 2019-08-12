context("Mortality")

library(raster)


# Create the samc object
samc_obj <- samc(res, abs, fid, tr_fun = function(x) 1/mean(x), override = TRUE)

# Extract Q
Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
Q <- as.matrix(Q)

# Extract R
R <- diag(nrow(Q))
diag(R) <- samc_obj@p[-nrow(samc_obj@p), ncol(samc_obj@p)]

# Create an indentity matrix
I <- diag(nrow(Q))

# Prepare the occupancy data
occ_ras <- raster(occ)
pv <- as.vector(occ_ras)
pv <- pv[is.finite(pv)]


# Run the tests
test_that("Testing mortality(samc, time)", {

  r1 <- mortality(samc_obj, time = time)

  r2 <- diag(nrow(Q))

  Qt <- diag(nrow(Q))

  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    r2 <- r2 + Qt
  }

  r2 <- r2 %*% R

  # Verify
  expect_equal(dim(r1), dim(r2))
  expect_equal(as.vector(r1), as.vector(r2))
})

test_that("Testing mortality(samc, origin, time)", {

  r1 <- mortality(samc_obj, origin = row, time = time)

  r2 <- diag(nrow(Q))

  Qt <- diag(nrow(Q))

  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    r2 <- r2 + Qt
  }

  r2 <- r2 %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, ]))
})

test_that("Testing mortality(samc, dest, time)", {

  r1 <- mortality(samc_obj, dest = col, time = time)

  r2 <- diag(nrow(Q))

  Qt <- diag(nrow(Q))

  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    r2 <- r2 + Qt
  }

  r2 <- r2 %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[, col]))
})

test_that("Testing mortality(samc, dest, time)", {

  r1 <- mortality(samc_obj, origin = row, dest = col, time = time)

  r2 <- diag(nrow(Q))

  Qt <- diag(nrow(Q))

  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    r2 <- r2 + Qt
  }

  r2 <- r2 %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, col]))
})

test_that("Testing mortality(samc, occ, time)", {

  r1 <- mortality(samc_obj, occ = occ, time = time)

  r2 <- I

  Qt <- diag(nrow(Q))

  for (i in 1:(time - 1)) {
    Qt <- Qt %*% Q
    r2 <- r2 + Qt
  }

  r2 <- pv %*% r2 %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})

test_that("Testing mortality(samc)", {

  r1 <- mortality(samc_obj)

  r2 <- solve(I - Q) %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})


test_that("Testing mortality(samc, origin)", {

  r1 <- mortality(samc_obj, origin = row)

  r2 <- solve(I - Q) %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, ]))
})

test_that("Testing mortality(samc, dest)", {

  r1 <- mortality(samc_obj, dest = col)

  r2 <- solve(I - Q) %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[, col]))
})

test_that("Testing mortality(samc, origin, dest)", {

  r1 <- mortality(samc_obj, origin = row, dest = col)

  r2 <- solve(I - Q) %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, col]))
})

test_that("Testing mortality(samc, occ)", {

  r1 <- mortality(samc_obj, occ = occ)

  r2 <- pv %*% solve(I - Q) %*% R

  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})
