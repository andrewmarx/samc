context("Dispersal")

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


#Run the tests
test_that("Testing dispersal(samc, dest, time)", {

  r1 <- dispersal(samc_obj, dest = col, time = time)

  qj <- Q[-col, col]

  Qj <- Q[-col,-col]

  Qji <- diag(nrow(Qj))
  r2 <- Qji

  for (i in 1:(time - 1)) {
    Qji <- Qji %*% Qj
    r2 <- r2 + Qji
  }

  r2 <- r2 %*% qj


  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})

test_that("Testing dispersal(samc, dest, time_vec)", {

  r1 <- dispersal(samc_obj, dest = col, time = time_vec)

  qj <- Q[-col, col]

  Qj <- Q[-col,-col]
  for (i in 1:length(time_vec)) {
    Qji <- diag(nrow(Qj))
    r2 <- Qji

    for (j in 1:(time_vec[i] - 1)) {
      Qji <- Qji %*% Qj
      r2 <- r2 + Qji
    }

    r2 <- r2 %*% qj

    # Verify
    expect_equal((r1[[i]]), as.vector(r2))
  }
})

test_that("Testing dispersal(samc, occ, dest, time)", {

  r1 <- dispersal(samc_obj, occ = occ, dest = col, time = time)

  qj <- Q[-col, col]

  Qj <- Q[-col,-col]

  Qji <- diag(nrow(Qj))
  r2 <- Qji

  for (i in 1:(time - 1)) {
    Qji <- Qji %*% Qj
    r2 <- r2 + Qji
  }

  r2 <- pv[-col] %*% (r2 %*% qj)

  # Verify
  expect_equal(r1, as.numeric(r2))
})

test_that("Testing dispersal(samc, occ, dest, time_vec)", {

  r1 <- dispersal(samc_obj, occ = occ, dest = col, time = time_vec)

  qj <- Q[-col, col]

  for (i in 1:length(time_vec)) {
    Qj <- Q[-col,-col]

    Qji <- diag(nrow(Qj))
    r2 <- Qji

    for (j in 1:(time_vec[i] - 1)) {
      Qji <- Qji %*% Qj
      r2 <- r2 + Qji
    }

    r2 <- pv[-col] %*% (r2 %*% qj)

    # Verify
    expect_equal(r1[[i]], as.numeric(r2))
  }
})

test_that("Testing dispersal(samc)", {

  r1 <- dispersal(samc_obj)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- (f - I) %*% fdg

  # Verify
  expect_equal(dim(r1), dim(r2))
  expect_equal(as.vector(r1), as.vector(r2))
})

# TODO Remove the skip once dispersal(samc, origin) is implemented
test_that("Testing dispersal(samc, origin)", {

  skip("dispersal(samc, origin) is not implemented")

  r1 <- dispersal(samc_obj, origin = row)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- (f - I) %*% fdg

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, ]))
})

test_that("Testing dispersal(samc, dest)", {

  r1 <- dispersal(samc_obj, dest = col)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- (f - I) %*% fdg

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[, col]))
})

test_that("Testing dispersal(samc, origin, dest)", {

  r1 <- dispersal(samc_obj, origin = row, dest = col)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- (f - I) %*% fdg

  # Verify
  expect_equal(as.vector(r1), as.vector(r2[row, col]))
})

test_that("Testing dispersal(samc, occ)", {

  r1 <- dispersal(samc_obj, occ = occ)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- pv %*% (f - I) %*% fdg

  # Verify
  expect_equal(as.vector(r1), as.vector(r2))
})


test_that("Testing dispersal(samc, occ, dest)", {

  r1 <- dispersal(samc_obj, occ = occ, dest = col)

  f <- solve(I - Q)

  fdg <- I
  diag(fdg) <- 1/diag(f)

  r2 <- pv %*% (f - I) %*% fdg

  # Verify
  expect_equal(r1, as.vector(r2)[col])
})
