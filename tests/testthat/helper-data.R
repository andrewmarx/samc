# testthat sources files beginning with helper- before running unit tests.
# This helper file contains data for each of the unit tests so that it is not
# repeated across test files.

baselist <- list(
  res = matrix(c( 1,  1,  2,  2,  3,  3,  2,  2,  1,
                  1,  2,  2,  3, 10, 10,  3,  2,  2,
                  1,  2,  2,  3, 10, 10,  3,  2,  2,
                  1,  1,  1,  3, 10, 10,  3,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1, 10, 10, 10, 10, 10, 10,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1),
                nrow = 9),
  abs = matrix(c( 1,  2,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  2,  2,  1,  1,  0,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  1,  5,  1,  1,  1,  1,  1,  0,  0,
                  1,  1,  1,  1,  1,  1,  1,  1,  0,
                  1,  1,  1,  2,  1,  1,  1,  7,  7,
                  3,  1,  1,  1,  1,  1,  1,  1,  1),
                nrow = 9) / 1000,
  fid = matrix(c(10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10,  1,  1, 10, 10, 10,
                 10, 10, 10, 10,  1,  1, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 50, 50, 50,
                 10, 10, 10, 10, 10, 10, 50, 50, 50,
                 10, 10, 10, 10, 10, 10, 10, 10, 10,
                 10, 10, 10, 10, 10, 10, 10, 10, 10),
                nrow = 9) / 100,
  init = matrix(c( 0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  1,  1,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  0),
                nrow = 9)
)


# Raster masks for testing different scenarios
masklist <- list(
  mask1 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask2 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1, NA, NA, NA, NA,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask3 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9),
  mask4 = matrix(c( 1,  1,  1,  1,  1,  1,  1,  1,  1,
                    1,  1,  1,  1, NA, NA, NA,  1,  1,
                    1,  1,  1,  1, NA,  1, NA,  1,  1,
                    1,  1,  1,  1, NA, NA, NA,  1,  1,
                    1, NA, NA, NA,  1,  1,  1,  1,  1,
                    1, NA,  1, NA,  1,  1,  1,  1,  1,
                    1, NA, NA, NA,  1,  1,  1,  1,  1,
                    1,  1,  1,  1,  1,  1,  1,  1,  1),
                  nrow = 9)
)

# Seed so the random absorption states (and any failures that depend on them)
# are reproducible from run to run.
set.seed(1)
p1 <- runif(72, max = 0.4)
p2 <- 1 - p1

testlist <- list()
for(i in seq_along(masklist)) {
  testlist[[i]] <- lapply(baselist, function(x) {masklist[[i]] * x})

  testlist[[i]]$length <- sum(!is.na(testlist[[i]]$res))

  testlist[[i]]$samc <- samc(testlist[[i]]$res,
                             testlist[[i]]$abs,
                             testlist[[i]]$fid,
                             model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))


  testlist[[i]]$samc@names = as.character(seq_along(testlist[[i]]$samc@data@t_abs))

  testlist[[i]]$samc$abs_states <- list(testlist[[i]]$abs * p1, testlist[[i]]$abs * p2)

  testlist[[i]]$id <- i
}

# Asymmetric versions
n <- length(testlist)
for(i in (n + 1):(n + length(masklist))) {
  testlist[[i]] <- lapply(baselist, function(x) {masklist[[i - n]] * x})

  testlist[[i]]$length <- sum(!is.na(testlist[[i]]$res))

  testlist[[i]]$samc <- samc(testlist[[i]]$res,
                             testlist[[i]]$abs,
                             testlist[[i]]$fid,
                             model = list(fun = function(x) 1/(mean(x) + x[1]), dir = 4, sym = FALSE))

  testlist[[i]]$samc@names = as.character(seq_along(testlist[[i]]$samc@data@t_abs))

  testlist[[i]]$samc$abs_states <- list(testlist[[i]]$abs * p1, testlist[[i]]$abs * p2)
  testlist[[i]]$id <- i
}


time = 100
time_vec = c(3, 5, 7, 11, 13)
row_vec = c(7, 34, 5, 5)
col_vec = c(13, 13, 5, 19)


#
# Shared reference-implementation helpers used across the metric tests. These
# build the "ground truth" results from first principles so they stay
# independent of the package's own implementation.
#

# M^t
mat_power <- function(M, t) {
  res <- M
  for (i in seq_len(t - 1)) {
    res <- res %*% M
  }
  res
}

# I + M + M^2 + ... + M^(t-1)
sum_powers <- function(M, t) {
  acc <- diag(nrow(M))
  term <- acc
  for (i in seq_len(t - 1)) {
    term <- term %*% M
    acc <- acc + term
  }
  acc
}

# Occupancy vector in the package's internal (row-major) cell ordering. The
# matrix is routed through raster() to match that ordering, then NA cells are
# dropped.
as_pv <- function(init) {
  v <- as.vector(raster::raster(init))
  v[is.finite(v)]
}

# Dense transient matrix Q for a samc object.
ref_q <- function(samc) {
  as.matrix(samc$q_matrix)
}

# Diagonal total-absorption matrix R = diag(t_abs).
ref_r <- function(samc) {
  R <- diag(length(samc@data@t_abs))
  diag(R) <- samc@data@t_abs
  R
}

# List of absorption matrices: the total R followed by one diagonal matrix per
# absorption state (the columns of c_abs). Used by the mortality tests.
ref_r_list <- function(samc) {
  per_state <- lapply(split(samc@data@c_abs, col(samc@data@c_abs)),
                      function(x) {
                        mat <- diag(length(x))
                        diag(mat) <- x
                        mat
                      })
  c(list(total = ref_r(samc)), per_state)
}

# Fundamental matrix F = (I - Q)^-1.
ref_f <- function(samc) {
  Q <- ref_q(samc)
  solve(diag(nrow(Q)) - Q)
}
