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

test_that("disconnected regions give a warning and connected inputs do not", {

  rw_model = list(fun = "1/mean(x)", dir = 8, sym = TRUE)

  res = matrix(1, 5, 5)
  abs = matrix(0.01, 5, 5)
  expect_no_warning(samc(res, abs, model = rw_model))

  # A column of NA cells splits the map into two regions
  res[, 3] = NA
  abs[, 3] = NA
  expect_warning(samc(res, abs, model = rw_model),
                 "Input contains disconnected regions")
})

test_that("P matrix construction gives a message and prints nothing", {

  p = matrix(c(0, 0.5, 0,
               0.5, 0, 0,
               0.5, 0.5, 1), 3, 3)

  expect_message(samc(p), "Some checks for manually created P matrices are still missing")
  expect_output(suppressMessages(samc(p)), NA)
})

test_that("P matrix construction rejects the conv method", {

  p = matrix(c(0, 0.5, 0,
               0.5, 0, 0,
               0.5, 0.5, 1), 3, 3)
  p_sparse = methods::as(methods::as(Matrix::Matrix(p, sparse = TRUE), "CsparseMatrix"), "generalMatrix")

  expect_error(samc(p, options = list(method = "conv")),
               "conv method is not supported for samc objects created from a P matrix")
  expect_error(samc(p_sparse, options = list(method = "conv")),
               "conv method is not supported for samc objects created from a P matrix")
})

test_that("non-square cells only give a warning when there is no CRS", {

  rw_model = list(fun = "1/mean(x)", dir = 8, sym = TRUE)
  res = matrix(seq(1, 2, length.out = 100), 10, 10)
  abs = matrix(0.01, 10, 10)
  ext = terra::ext(0, 20, 0, 10)

  res_crs = terra::rast(res, extent = ext, crs = "EPSG:32617")
  abs_crs = terra::rast(abs, extent = ext, crs = "EPSG:32617")
  expect_no_warning(samc(res_crs, abs_crs, model = rw_model))

  res_nocrs = terra::rast(res, extent = ext, crs = "")
  abs_nocrs = terra::rast(abs, extent = ext, crs = "")
  expect_warning(samc(res_nocrs, abs_nocrs, model = rw_model),
                 "Raster cells are not square")
})

test_that("negative data values give an error", {

  res = matrix(1, 5, 5)
  res[2, 2] = -1
  abs = matrix(0.01, 5, 5)

  expect_error(samc(res, abs, model = list(fun = "1/mean(x)", dir = 8, sym = TRUE)),
               "The data must not have negative values")
})

test_that("an invalid dir gives an error naming the argument", {
  res = testlist[[1]]$res
  abs = testlist[[1]]$abs
  msg = "`dir` must be set to either 4 or 8"

  models = list(
    function(dir) list(fun = "1/mean(x)", dir = dir, sym = TRUE),
    function(dir) list(name = "CRW", fun = "1/mean(x)", dir = dir, sym = TRUE, dist = "vonMises", kappa = 1),
    function(dir) list(name = "SSF", fun = "x[2]", dir = dir, sym = FALSE, ssc = 1)
  )

  for (model in models) {
    expect_error(samc(res, abs, model = model(5)), msg)
    expect_error(samc(res, abs, model = model(0)), msg)
    expect_error(samc(res, abs, model = model(c(4, 8))), msg)
    expect_error(samc(res, abs, model = model(NA)), msg)
    expect_error(samc(res, abs, model = model(NA_real_)), msg)
    expect_error(samc(res, abs, model = model("8")), msg)
    expect_error(samc(res, abs, model = model(numeric(0))), msg)
    expect_error(samc(res, abs, model = model(NULL)), msg)
    expect_error(samc(res, abs, model = model(list(8))), msg)
    expect_error(samc(res, abs, model = model(TRUE)), msg)
    expect_error(samc(res, abs, model = model(matrix(8))), msg)
  }
})

test_that("missing, unknown, duplicate, and unnamed model arguments are listed in a single message", {
  res = testlist[[1]]$res
  abs = testlist[[1]]$abs

  expect_error(samc(res, abs, model = list(dir = 8, sym = TRUE)), "Missing argument in model: fun$")
  expect_error(samc(res, abs, model = list(sym = TRUE)), "Missing argument in model: fun, dir$")
  expect_error(samc(res, abs, model = list(name = "SSF", fun = "x[2]")), "Missing argument in model: dir, sym, ssc$")
  expect_error(samc(res, abs, model = list(name = "CRW", fun = "1/mean(x)", dir = 8, sym = TRUE)),
               "Missing argument in model: dist, kappa$")

  expect_error(samc(res, abs, model = list(fun = "1/mean(x)", dir = 8, sym = TRUE, foo = 1)),
               "Unknown argument in model: foo$")
  expect_error(samc(res, abs, model = list(foo = 1, fun = "1/mean(x)", dir = 8, sym = TRUE, bar = 2)),
               "Unknown argument in model: foo, bar$")

  expect_error(samc(res, abs, model = list(fun = "1/mean(x)", dir = 8, dir = 4, sym = TRUE)),
               "Duplicate argument in model: dir$")
  expect_error(samc(res, abs, model = list(fun = "1/mean(x)", dir = 8, dir = 4, dir = 4, sym = TRUE, sym = TRUE)),
               "Duplicate argument in model: dir, sym$")

  expect_error(samc(res, abs, model = list("1/mean(x)", 8, TRUE)), "model: all elements must be named")
  expect_error(samc(res, abs, model = list(fun = "1/mean(x)", 8, sym = TRUE)), "model: all elements must be named")
})

test_that("the conv method only accepts the 1/mean(x) named function", {
  res = testlist[[1]]$res
  abs = testlist[[1]]$abs
  msg = "Convolution currently only supports the '1/mean\\(x\\)' named function"
  conv = list(method = "conv")

  expect_error(samc(res, abs, model = list(fun = function(x) 1 / mean(x), dir = 8, sym = TRUE), options = conv), msg)
  expect_no_error(suppressWarnings(samc(res, abs, model = list(fun = "1/mean(x)", dir = 8, sym = TRUE), options = conv)))
})

test_that("an invalid fun gives an error naming the argument and the supported names", {
  res = testlist[[1]]$res
  abs = testlist[[1]]$abs

  rw = function(fun) list(fun = fun, dir = 8, sym = TRUE)
  crw = function(fun) list(name = "CRW", fun = fun, dir = 8, sym = TRUE, dist = "vonMises", kappa = 1)
  ssf = function(fun) list(name = "SSF", fun = fun, dir = 8, sym = FALSE, ssc = 1)

  expect_error(samc(res, abs, model = rw("x[2]")),
               "`fun` = \"x\\[2\\]\" is not a named function supported by the RW model. Supported names: \"1/mean\\(x\\)\"$")
  expect_error(samc(res, abs, model = crw("x[2]")),
               "`fun` = \"x\\[2\\]\" is not a named function supported by the CRW model. Supported names: \"1/mean\\(x\\)\"$")
  expect_error(samc(res, abs, model = ssf("1/mean(x)")),
               "`fun` = \"1/mean\\(x\\)\" is not a named function supported by the SSF model. Supported names: \"x\\[2\\]\"$")
  expect_error(samc(res, abs, model = rw("foo")), "`fun` = \"foo\" is not a named function supported by the RW model")

  # Reported before the method check
  expect_error(samc(res, abs, model = rw("x[2]"), options = list(method = "conv")), "is not a named function supported by the RW model")

  msg = "`fun` must be a function or the name of a supported named function"
  for (model in list(rw, crw, ssf)) {
    expect_error(samc(res, abs, model = model(1)), msg)
    expect_error(samc(res, abs, model = model(NULL)), msg)
    expect_error(samc(res, abs, model = model(NA_character_)), msg)
    expect_error(samc(res, abs, model = model(c("1/mean(x)", "x[2]"))), msg)
    expect_error(samc(res, abs, model = model(list("1/mean(x)"))), msg)
  }
})
