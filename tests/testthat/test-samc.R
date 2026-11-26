# The expected row sums: one entry per non-NA cell plus the extra absorbing row
# of the p matrix, each summing to 1.
expected_rowsums <- rep(1, sum(!is.na(as.vector(testlist[[1]]$res))) + 1)

test_that("samc p matrix row sums equal 1 without fidelity data", {

  # Create the samc object (no fidelity) and get the row sums of the p matrix
  samc_obj <- samc(testlist[[1]]$res, testlist[[1]]$abs, model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))
  rs <- Matrix::rowSums(samc_obj$p_matrix)

  # Verify equality
  expect_equal(rs, expected_rowsums, check.names = FALSE)
})

test_that("samc p matrix row sums equal 1 with fidelity data", {

  # The shared fixture is built with fidelity data, so reuse it here
  rs <- Matrix::rowSums(testlist[[1]]$samc$p_matrix)

  # Verify equality
  expect_equal(rs, expected_rowsums, check.names = FALSE)
})

test_that("samc objects built from P matrices do not share a cache", {

  # Two absorbing chains with the same structure but different transition
  # probabilities. Solving one must not affect the results of the other.
  make_p <- function(q) {
    p <- matrix(0, 3, 3)
    p[1, ] <- c(0, q, 1 - q)
    p[2, ] <- c(q, 0, 1 - q)
    p[3, 3] <- 1
    methods::as(methods::as(Matrix::Matrix(p, sparse = TRUE), "CsparseMatrix"), "generalMatrix")
  }

  samc_a = suppressMessages(samc(make_p(0.5)))
  samc_b = suppressMessages(samc(make_p(0.9)))

  expect_false(identical(samc_a@.cache, samc_b@.cache))

  # Solve a first so that a shared cache would hold its factorization
  vis_a <- visitation(samc_a, origin = 1)
  vis_b <- visitation(samc_b, origin = 1)

  expect_equal(vis_a, as.vector(solve(diag(2) - as.matrix(samc_a$q_matrix))[1, ]), tolerance = 1e-10)
  expect_equal(vis_b, as.vector(solve(diag(2) - as.matrix(samc_b$q_matrix))[1, ]), tolerance = 1e-10)
})

test_that("SpatRasters without a CRS work with the convolution and CRW paths", {

  # `terra::is.lonlat()` returns NA (with a warning) when a raster has no CRS.
  # That must not break the lonlat checks in samc().
  res = matrix(seq(1, 2, length.out = 100), 10, 10)
  abs = matrix(0.01, 10, 10)
  res_r = terra::rast(res)
  abs_r = terra::rast(abs)
  expect_true(is.na(suppressWarnings(terra::is.lonlat(res_r))))

  rw_model = list(fun = "1/mean(x)", dir = 8, sym = TRUE)
  conv_opts = list(method = "conv", threads = 1)

  conv_r = suppressWarnings(samc(res_r, abs_r, model = rw_model, options = conv_opts))
  conv_m = samc(res, abs, model = rw_model, options = conv_opts)
  expect_equal(mortality(conv_r, init = res_r, time = 5),
               mortality(conv_m, init = res, time = 5),
               tolerance = 1e-12)

  crw_model = list(name = "CRW", fun = "1/mean(x)", dir = 8, sym = TRUE, dist = "vonMises", kappa = 1)
  expect_no_error(suppressWarnings(samc(res_r, abs_r, model = crw_model)))
})

test_that("lonlat warnings still fire for the convolution and CRW paths", {

  res_r = terra::rast(matrix(seq(1, 2, length.out = 100), 10, 10),
                      extent = terra::ext(0, 10, 0, 10), crs = "EPSG:4326")
  abs_r = terra::rast(matrix(0.01, 10, 10),
                      extent = terra::ext(0, 10, 0, 10), crs = "EPSG:4326")

  collect_warnings = function(expr) {
    msgs = character(0)
    withCallingHandlers(expr, warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    msgs
  }

  conv_warnings = collect_warnings(samc(res_r, abs_r, model = list(fun = "1/mean(x)", dir = 8, sym = TRUE),
                                        options = list(method = "conv", threads = 1)))
  expect_true(any(grepl("Convolution method currently does not properly correct directions for lonlat", conv_warnings)))

  crw_warnings = collect_warnings(samc(res_r, abs_r, model = list(name = "CRW", fun = "1/mean(x)", dir = 8, sym = TRUE,
                                                                  dist = "vonMises", kappa = 1)))
  expect_true(any(grepl("CRW does not properly adjust turning angles for lonlat", crw_warnings)))
})
