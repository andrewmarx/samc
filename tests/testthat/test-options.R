# Tests for validation of the `options` and `model` arguments of samc()

opt_res = samc::example_split_corridor$res
opt_abs = samc::example_split_corridor$abs
opt_model = list(fun = "1/mean(x)", dir = 8, sym = TRUE)

test_that("method 'iter' is accepted", {
  samc_obj = suppressWarnings(samc(opt_res, opt_abs, model = opt_model, options = list(method = "iter")))
  expect_equal(samc_obj$solver, "iter")

  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = "foo")),
               "method must be 'direct', 'iter', or 'conv'")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = c("direct", "iter"))),
               "method must be 'direct', 'iter', or 'conv'")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = NA_character_)),
               "method must be 'direct', 'iter', or 'conv'")
})

test_that("unknown option names are reported by name", {
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(foo = 1)),
               "Invalid option names: foo$")
  expect_error(samc(opt_res, opt_abs, model = opt_model,
                    options = list(method = "direct", threads = 1, override = FALSE, precision = "double", foo = 1)),
               "Invalid option names: foo$")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(foo = 1, method = "direct", bar = 2)),
               "Invalid option names: foo bar$")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list("direct")),
               "all elements must be named")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = "direct", 1)),
               "all elements must be named")
})

test_that("invalid precision is reported as precision", {
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(precision = "half")),
               "precision must be 'single' or 'double'")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(precision = NA_character_)),
               "precision must be 'single' or 'double'")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(precision = c("single", "double"))),
               "precision must be 'single' or 'double'")
})

test_that("single precision is rejected for methods other than conv", {
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = "direct", precision = "single")),
               "method: 'direct' only supports precision: 'double'")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(method = "iter", precision = "single")),
               "method: 'iter' only supports precision: 'double'")
})

test_that("override must be a single non-NA logical", {
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(override = NA)),
               "override must be TRUE or FALSE")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(override = c(TRUE, FALSE))),
               "override must be TRUE or FALSE")
  expect_error(samc(opt_res, opt_abs, model = opt_model, options = list(override = "TRUE")),
               "override must be TRUE or FALSE")
})

test_that("sym must be a single non-NA logical", {
  expect_error(samc(opt_res, opt_abs, model = list(fun = "1/mean(x)", dir = 8, sym = NA)),
               "`sym` must be set to either TRUE or FALSE")
  expect_error(samc(opt_res, opt_abs, model = list(fun = "1/mean(x)", dir = 8, sym = c(TRUE, FALSE))),
               "`sym` must be set to either TRUE or FALSE")
})
