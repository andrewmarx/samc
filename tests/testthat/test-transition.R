# Tests for the transition values of user defined transition functions

test_that("edges from or to NA cells or off the map have NA transition values for every function", {
    m = matrix(c(1, 2, 3, 4, 5, NA, 7, 8, 9, 10, 11, 12), nrow = 3, byrow = TRUE)
    abs = m * 0 + 0.05
    vals = as.vector(t(m)) # cell order
    r = terra::rast(m, crs = "local")

    const = function(x) 1

    offsets = list(
        "4" = data.frame(dr = c(-1, 0, 0, 1), dc = c(0, -1, 1, 0)),
        "8" = data.frame(dr = c(-1, -1, -1, 0, 0, 1, 1, 1), dc = c(-1, 0, 1, -1, 1, -1, 0, 1))
    )

    for (dir in c(4, 8)) {
        # The same NA pattern as the named function, for RW/CRW and SSF values
        named_na = is.na(samc:::.tr_vals(r, "1/mean(x)", dir))
        expect_identical(is.na(samc:::.tr_vals(r, const, dir)), named_na)
        expect_identical(is.na(samc:::.tr_vals(r, const, dir, ssc = 1.5)), named_na)

        # Explicit Q matrix: the value of each edge between two non-NA cells is 1 / distance, normalized over those
        # edges only
        cells = which(!is.na(vals))
        q = matrix(0, length(cells), length(cells))
        for (s in seq_along(cells)) {
            row = (cells[s] - 1) %/% ncol(m) + 1
            col = (cells[s] - 1) %% ncol(m) + 1
            t = numeric(length(cells))
            for (d in seq_len(dir)) {
                dr = offsets[[as.character(dir)]]$dr[d]
                dc = offsets[[as.character(dir)]]$dc[d]
                r2 = row + dr
                c2 = col + dc
                if (r2 < 1 || r2 > nrow(m) || c2 < 1 || c2 > ncol(m)) next
                to = (r2 - 1) * ncol(m) + c2
                if (is.na(vals[to])) next
                t[match(to, cells)] = 1 / sqrt(dr^2 + dc^2)
            }
            q[s, ] = (1 - 0.05) * t / sum(t)
        }

        samc_obj = suppressWarnings(samc(m, abs, model = list(fun = const, dir = dir, sym = TRUE)))
        expect_equal(unname(as.matrix(samc_obj$q_matrix)), q, tolerance = 1e-12)
        expect_equal(unname(rowSums(as.matrix(samc_obj$q_matrix))), rep(0.95, length(cells)), tolerance = 1e-12)
    }

    # Global lonlat raster with NA cells: wrapped edges are kept, others follow the named function
    g = terra::rast(nrows = 4, ncols = 6, xmin = -180, xmax = 180, ymin = -60, ymax = 60, crs = "EPSG:4326")
    terra::values(g) = c(1:7, NA, 9:24)
    for (dir in c(4, 8)) {
        expect_identical(is.na(samc:::.tr_vals(g, const, dir)), is.na(samc:::.tr_vals(g, "1/mean(x)", dir)))
    }
})

test_that("functions that ignore NA give the same results as NA aware functions", {
    na_aware = function(x) if (anyNA(x)) NA else 1

    res = testlist[[2]]$res
    abs = testlist[[2]]$abs

    rw = function(f) list(fun = f, dir = 8, sym = TRUE)
    a = suppressWarnings(samc(res, abs, model = rw(function(x) 1)))
    b = suppressWarnings(samc(res, abs, model = rw(na_aware)))
    expect_identical(a@data@f, b@data@f)

    ssf = function(f) list(name = "SSF", fun = f, dir = 8, sym = FALSE, ssc = 1.5)
    a = suppressWarnings(samc(res, abs, model = ssf(function(x) 1)))
    b = suppressWarnings(samc(res, abs, model = ssf(na_aware)))
    expect_identical(a@data@f, b@data@f)

    # The CRW model does not support NA cells, so this only covers the edges of the map
    res = testlist[[1]]$res
    abs = testlist[[1]]$abs
    crw = function(f) list(name = "CRW", fun = f, dir = 8, sym = TRUE, dist = "vonMises", kappa = 1)
    a = suppressWarnings(samc(res, abs, model = crw(function(x) 1)))
    b = suppressWarnings(samc(res, abs, model = crw(na_aware)))
    expect_identical(a@data@f, b@data@f)
    expect_identical(a@crw_map, b@crw_map)
    expect_identical(a@prob_mat, b@prob_mat)
})
