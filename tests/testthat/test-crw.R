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

# Long-term dispersal to a dest cell. The reference makes all states of the dest
# cell absorbing and solves for the probability of ever entering them. For the
# dest cell itself the value is the probability of returning to it.
crw_disp_r = terra::rast(matrix(c(3, 7, 1, 9, 4, 2, 8, 5, 6, 1, 4, 9, 2, 7, 3, 5, 1, 6, 8, 2, 9, 3, 7, 4, 6), 5, 5),
                         crs = "local")
crw_disp_fid = terra::rast(matrix(0.1, 5, 5), crs = "local")

# Per state (row of crw_map)
crw_disp_ref_states = function(samc_obj, dest) {
  q = as.matrix(samc_obj$q_matrix)
  S = samc_obj@crw_map[, 1] == dest

  h = numeric(nrow(q))
  h[!S] = solve(diag(sum(!S)) - q[!S, !S], rowSums(q[!S, S, drop = FALSE]))
  h[S] = rowSums(q[S, S, drop = FALSE]) + q[S, !S] %*% h[!S]
  h
}

# Per cell, for a walker whose heading in each cell follows samc_obj@prob_mat
crw_disp_ref = function(samc_obj, dest) {
  pv = samc_obj@prob_mat
  pv = pv[!is.na(pv)]

  as.vector(tapply(pv * crw_disp_ref_states(samc_obj, dest), samc_obj@crw_map[, 1], sum))
}

for (dir in c(4, 8)) {
  for (kappa in c(0, 2, 10)) {
    test_that(paste0("CRW long-term dispersal(samc, dest) is exact (dir = ", dir, ", kappa = ", kappa, ")"), {
      samc_obj = suppressWarnings(samc(crw_disp_r, crw_abs, crw_disp_fid,
        model = list(name = "CRW", fun = "1/mean(x)", dir = dir, sym = TRUE, dist = "vonMises", kappa = kappa)))

      for (dest in c(1, 8, 13)) {
        ref = crw_disp_ref(samc_obj, dest)

        expect_equal(as.vector(dispersal(samc_obj, dest = dest)), ref, tolerance = 1e-10)

        # Off the dest cell, the long-term result is the limit of the time-step result
        expect_equal(dispersal(samc_obj, dest = dest, time = 20000)[-dest], ref[-dest], tolerance = 1e-10)

        expect_equal(dispersal(samc_obj, origin = 7, dest = dest), ref[7], tolerance = 1e-10)
      }

      # init is per state, with a heading mix that differs from prob_mat
      set.seed(1)
      init = runif(nrow(samc_obj@crw_map))
      expect_equal(dispersal(samc_obj, init = init, dest = 13), sum(init * crw_disp_ref_states(samc_obj, 13)),
                   tolerance = 1e-10)

      samc_obj$solver = "iter"
      expect_equal(as.vector(dispersal(samc_obj, dest = 13)), crw_disp_ref(samc_obj, 13), tolerance = 1e-6)
    })
  }
}

test_that("CRW long-term dispersal() without dest is disabled", {
  samc_obj = suppressWarnings(samc(crw_disp_r, crw_abs, crw_disp_fid,
    model = list(name = "CRW", fun = "1/mean(x)", dir = 4, sym = TRUE, dist = "vonMises", kappa = 2)))

  expect_error(dispersal(samc_obj, origin = 7), "not currently supported for CRW")
  expect_error(dispersal(samc_obj, init = rep(1, nrow(samc_obj@crw_map))), "not currently supported for CRW")
})
