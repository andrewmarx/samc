context("Conditional Passage Time")


# Create the samc object
samc_obj <- samc(res, abs, fid, tr_fun = function(x) 1/mean(x), override = TRUE)

# Create a version from P matrix
samc_p <- samc(p_mat = samc_obj@p)


# Calculate the results based on De Sanctis and de Koning 2018
Q <- samc_obj@p[-nrow(samc_obj@p), -ncol(samc_obj@p)]
Q <- as.matrix(Q)

qj <- Q[-col, col]
Qj <- Q[-col, -col]

I <- diag(nrow(Qj))

r <- samc_obj@p[-nrow(samc_obj@p), ncol(samc_obj@p)]
r <- r[-col]

R <- cbind(r, qj)

f <- solve(I - Qj)

b <- as.matrix(f %*% R)
bdg <- Matrix::sparseMatrix(i = 1:nrow(b),
                            j = 1:nrow(b),
                            x = b[, 2],
                            index1 = TRUE)

result <- solve(bdg) %*% f %*% bdg %*% rep(1, nrow(bdg))
result <- as.numeric(result)

# Run the tests
test_that("Testing cond_passage(samc, dest)", {

  expect_error(cond_passage(samc_obj, dest = col))

  r1 <- cond_passage(samc_p, dest = col)

  # Verify
  expect_equal(dim(r1), dim(result))
  expect_equal(as.vector(r1), as.vector(result))
})

test_that("Testing cond_passage(samc, origin, dest)", {

  expect_error(cond_passage(samc_obj, origin = row, dest = col))

  r1 <- cond_passage(samc_p, origin = row, dest = col)

  # Verify
  expect_equal(as.vector(r1), as.vector(result[row]))
})
