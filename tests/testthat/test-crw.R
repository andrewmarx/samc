# Tests for the CRW turning kernel. Heading indices in `crw_map[, 2]` follow
# the direction order of the model: dir = 4 is N, W, E, S and dir = 8 is
# NW, N, NE, W, E, SW, S, SE.

crw_r = terra::rast(matrix(1, 5, 5), crs = "local")
crw_abs = terra::rast(matrix(0.01, 5, 5), crs = "local")
crw_kappa = 2

# P(next heading | current heading) for the edges leaving the centre cell,
# restricted to the 4 compass headings (N, W, E, S) and row-normalised.
crw_turn_probs = function(dir) {
  samc_obj = suppressWarnings(samc(crw_r, crw_abs,
    model = list(name = "CRW", fun = "1/mean(x)", dir = dir, sym = TRUE, dist = "vonMises", kappa = crw_kappa)))

  compass = if (dir == 4) 1:4 else c(2, 4, 5, 7)
  crw_map = samc_obj@crw_map
  q = as.matrix(samc_obj$q_matrix)

  centre = 13
  probs = matrix(NA_real_, 4, 4, dimnames = list(c("N", "W", "E", "S"), c("N", "W", "E", "S")))
  for (i in 1:4) {
    row = q[crw_map[, 1] == centre & crw_map[, 2] == compass[i], ]
    for (j in 1:4) {
      probs[i, j] = sum(row[crw_map[, 2] == compass[j]])
    }
  }

  probs / rowSums(probs)
}

# Von Mises weights over the 4 compass headings: cos = 1 straight, 0 for a side
# turn, -1 for a reversal
crw_cos = matrix(c( 1,  0,  0, -1,
                    0,  1, -1,  0,
                    0, -1,  1,  0,
                   -1,  0,  0,  1), 4, 4, byrow = TRUE)
crw_expected = exp(crw_kappa * crw_cos)
crw_expected = crw_expected / rowSums(crw_expected)
dimnames(crw_expected) = list(c("N", "W", "E", "S"), c("N", "W", "E", "S"))

test_that("CRW turning kernel follows the von Mises distribution for dir = 4", {
  expect_equal(crw_turn_probs(4), crw_expected, tolerance = 1e-10)
})

test_that("CRW turning kernel follows the von Mises distribution for dir = 8", {
  expect_equal(crw_turn_probs(8), crw_expected, tolerance = 1e-10)
})

# P(next heading | current heading) for a walker heading N from cell 8 into the
# top-row cell 3, where heading N is unavailable. Returned over the headings
# that remain, row-normalised.
crw_edge_probs = function(dir, kappa) {
  samc_obj = suppressWarnings(samc(crw_r, crw_abs,
    model = list(name = "CRW", fun = "1/mean(x)", dir = dir, sym = TRUE, dist = "vonMises", kappa = kappa)))

  north = if (dir == 4) 1 else 2
  nxt = if (dir == 4) c(W = 2, E = 3, S = 4) else c(W = 4, E = 5, SW = 6, S = 7, SE = 8)
  crw_map = samc_obj@crw_map
  q = as.matrix(samc_obj$q_matrix)

  row = q[crw_map[, 1] == 8 & crw_map[, 2] == north, ]
  probs = sapply(nxt, function(h) sum(row[crw_map[, 1] == 3 & crw_map[, 2] == h]))

  probs / sum(probs)
}

crw_edge_expected = function(dir, kappa) {
  cos_turn = if (dir == 4) c(W = 0, E = 0, S = -1) else c(W = 0, E = 0, SW = -sqrt(0.5), S = -1, SE = -sqrt(0.5))
  # Diagonal steps are down-weighted by their length in the transition values
  step_wt = if (dir == 4) c(1, 1, 1) else c(1, 1, sqrt(0.5), 1, sqrt(0.5))
  w = step_wt * exp(kappa * cos_turn)
  w / sum(w)
}

for (dir in c(4, 8)) {
  for (kappa in c(2, 500, 1e5, 1e6)) {
    test_that(paste0("CRW kernel is valid at the map edge (dir = ", dir, ", kappa = ", kappa, ")"), {
      expect_equal(crw_edge_probs(dir, kappa), crw_edge_expected(dir, kappa), tolerance = 1e-10)
    })
  }
}

test_that("CRW with large kappa builds a valid P matrix", {
  for (dir in c(4, 8)) {
    for (kappa in c(500, 5000, 1e5, 1e6)) {
      samc_obj = suppressWarnings(samc(crw_r, crw_abs,
        model = list(name = "CRW", fun = "1/mean(x)", dir = dir, sym = TRUE, dist = "vonMises", kappa = kappa)))

      expect_false(anyNA(samc_obj@data@f@x))
      expect_equal(as.vector(survival(samc_obj)), rep(100, 25), tolerance = 1e-6)
    }
  }
})
