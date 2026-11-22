# Tests for global lonlat rasters, where terra::focalValues() wraps columns
# across the antimeridian. Edges built from those values must wrap within the
# same row.

global_r = terra::rast(nrows = 6, ncols = 12, xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = "EPSG:4326")
terra::values(global_r) = 1

global_abs = global_r
terra::values(global_abs) = 0.01

global_nc = terra::ncol(global_r)

# Row/col steps for each direction, in .tr_vals() order
global_steps = list(
  "4" = list(dr = c(-1, 0, 0, 1), dc = c(0, -1, 1, 0)),
  "8" = list(dr = c(-1, -1, -1, 0, 0, 1, 1, 1), dc = c(-1, 0, 1, -1, 1, -1, 0, 1))
)

test_that("neighbor offsets do not wrap on non-global rasters", {
  r = terra::rast(nrows = 4, ncols = 5, crs = "local")

  expect_equal(.nb_offsets(r, 4), matrix(c(-5, -1, 1, 5), 5, 4, byrow = TRUE))
  expect_equal(.nb_offsets(r, 8), matrix(c(-6, -5, -4, -1, 1, 4, 5, 6), 5, 8, byrow = TRUE))

  regional = terra::rast(nrows = 4, ncols = 5, xmin = -20, xmax = 20, ymin = 30, ymax = 60, crs = "EPSG:4326")
  expect_equal(.nb_offsets(regional, 8), matrix(c(-6, -5, -4, -1, 1, 4, 5, 6), 5, 8, byrow = TRUE))
})

for (dir in c(4, 8)) {
  steps = global_steps[[as.character(dir)]]

  test_that(paste("neighbor offsets wrap on global lonlat rasters for dir =", dir), {
    off = .nb_offsets(global_r, dir)
    row = 3

    for (col in c(1, 2, global_nc - 1, global_nc)) {
      cell = terra::cellFromRowCol(global_r, row, col)
      expected = terra::cellFromRowCol(global_r, row + steps$dr, (col + steps$dc - 1) %% global_nc + 1)
      expect_equal(cell + off[col, ], expected)
    }
  })

  test_that(paste("global lonlat samc wraps edges at the antimeridian for dir =", dir), {
    invisible(utils::capture.output(
      samc_obj <- samc(global_r, global_abs, model = list(fun = "1/mean(x)", dir = dir, sym = TRUE))
    ))
    q = as.matrix(samc_obj$q_matrix)

    # Every cell's edges point at its wrapped neighbors, including the top and bottom rows
    nr = terra::nrow(global_r)
    for (row in c(1, 3, nr)) {
      for (col in c(1, global_nc)) {
        cell = terra::cellFromRowCol(global_r, row, col)
        nb_rows = row + steps$dr
        keep = nb_rows >= 1 & nb_rows <= nr
        nbs = terra::cellFromRowCol(global_r, nb_rows[keep], (col + steps$dc[keep] - 1) %% global_nc + 1)

        expect_setequal(which(q[cell, ] > 0), nbs)
      }
    }

    expect_equal(unname(rowSums(q)), rep(0.99, nrow(q)), tolerance = 1e-12)
    expect_true(isSymmetric(1 * (q > 0)))

    # A uniform global raster is invariant to rotation in longitude, so results
    # from a seam cell must match those from an interior cell shifted by half
    # the columns
    half = global_nc / 2
    v_seam = visitation(samc_obj, origin = terra::cellFromRowCol(global_r, 3, 1))
    v_mid = visitation(samc_obj, origin = terra::cellFromRowCol(global_r, 3, 1 + half))

    m_seam = matrix(v_seam, ncol = global_nc, byrow = TRUE)
    m_mid = matrix(v_mid, ncol = global_nc, byrow = TRUE)

    expect_equal(m_seam, m_mid[, c((half + 1):global_nc, 1:half)], tolerance = 1e-10)
  })
}

test_that("global lonlat CRW builds an edge for every wrapped RW edge", {
  invisible(utils::capture.output(
    rw <- samc(global_r, global_abs, model = list(fun = "1/mean(x)", dir = 8, sym = TRUE))
  ))

  crw = suppressWarnings(samc(global_r, global_abs,
    model = list(name = "CRW", fun = "1/mean(x)", dir = 8, sym = TRUE, dist = "vonMises", kappa = 1)))

  expect_equal(nrow(crw$q_matrix), sum(rw$q_matrix > 0))
  expect_equal(unname(Matrix::rowSums(crw$q_matrix)), rep(0.99, nrow(crw$q_matrix)), tolerance = 1e-10)
})
