# Tests for the step length lookup used when building transition matrices.
# Columns of the lookup are 3x3 focal window positions, row-major from the
# top-left: 1 = NW, 2 = N, 3 = NE, 4 = W, 5 = center, 6 = E, 7 = SW, 8 = S, 9 = SE.

lookup_win = list("4" = c(2, 4, 6, 8), "8" = c(1:4, 6:9))
lookup_dr = rep(c(-1, 0, 1), each = 3)
lookup_dc = rep(c(-1, 0, 1), times = 3)

# Regional lonlat raster (does not span the globe, so no column wrapping)
lookup_ll = terra::rast(nrows = 6, ncols = 8, xmin = -20, xmax = 20, ymin = 30, ymax = 60, crs = "EPSG:4326")
terra::values(lookup_ll) = 1

test_that("planar step lengths are 1 and sqrt(2)", {
  r = terra::rast(nrows = 4, ncols = 5, crs = "local")

  m4 = .build_lookup_mat(r, 4)
  expect_equal(unique(as.vector(m4[, c(2, 4, 6, 8)])), 1)
  expect_true(all(is.na(m4[, c(1, 3, 5, 7, 9)])))

  m8 = .build_lookup_mat(r, 8)
  expect_equal(unique(as.vector(m8[, c(2, 4, 6, 8)])), 1)
  expect_equal(unique(as.vector(m8[, c(1, 3, 7, 9)])), sqrt(2))
})

for (dir in c(4, 8)) {
  test_that(paste("lonlat step lengths match cell-to-cell distances for dir =", dir), {
    expect_no_warning(m <- .build_lookup_mat(lookup_ll, dir))

    nr = terra::nrow(lookup_ll)
    col = 4 # interior column so every horizontal neighbor exists

    for (d in lookup_win[[as.character(dir)]]) {
      for (row in seq_len(nr)) {
        nb_row = row + lookup_dr[d]

        if (nb_row < 1 || nb_row > nr) {
          expect_true(is.na(m[row, d]))
        } else {
          from = terra::xyFromCell(lookup_ll, terra::cellFromRowCol(lookup_ll, row, col))
          to = terra::xyFromCell(lookup_ll, terra::cellFromRowCol(lookup_ll, nb_row, col + lookup_dc[d]))
          expect_equal(m[row, d], terra::distance(from, to, lonlat = TRUE, pairwise = TRUE), tolerance = 1e-8)
        }
      }
    }

    # Columns for unused directions are left empty
    expect_true(all(is.na(m[, setdiff(1:9, lookup_win[[as.character(dir)]])])))
  })

  test_that(paste("lonlat samc keeps edges in every direction for dir =", dir), {
    abs_r = lookup_ll
    terra::values(abs_r) = 0.01

    invisible(utils::capture.output(samc_obj <- samc(lookup_ll, abs_r, model = list(fun = "1/mean(x)", dir = dir, sym = TRUE))))
    q = samc_obj$q_matrix

    # Interior cell: every neighbor in the requested directions is reachable
    row = 3
    col = 4
    cell = terra::cellFromRowCol(lookup_ll, row, col)
    nbs = sapply(lookup_win[[as.character(dir)]], function(d) {
      terra::cellFromRowCol(lookup_ll, row + lookup_dr[d], col + lookup_dc[d])
    })

    expect_true(all(q[cell, nbs] > 0))
    expect_equal(sum(q[cell, ] > 0), dir)
  })
}
