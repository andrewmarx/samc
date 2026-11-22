# Tests for changing the solver of an existing samc object with `$solver<-`

test_that("solver can be switched between direct and iter", {
  samc_obj = testlist[[1]]$samc
  expect_equal(samc_obj$solver, "direct")

  vis_direct = visitation(samc_obj, origin = 1)

  samc_obj$solver = "iter"
  expect_equal(samc_obj$solver, "iter")
  expect_equal(visitation(samc_obj, origin = 1), vis_direct, tolerance = 1e-8)

  samc_obj$solver = "direct"
  expect_equal(samc_obj$solver, "direct")
  expect_equal(visitation(samc_obj, origin = 1), vis_direct, tolerance = 1e-12)
})

test_that("invalid solver values warn and leave the solver unchanged", {
  samc_obj = testlist[[1]]$samc

  for (value in list("conv", "foo", c("direct", "iter"), 1, NA_character_)) {
    expect_warning(samc_obj$solver <- value, "Valid inputs")
    expect_equal(samc_obj$solver, "direct")
  }
})

test_that("solver of a convolution samc object cannot be changed", {
  conv_obj = samc(testlist[[1]]$res,
                   pmax(testlist[[1]]$abs, 0.00001),
                   testlist[[1]]$fid,
                   model = list(fun = "1/mean(x)", dir = 8, sym = TRUE),
                   options = list(threads = 1, method = "conv", override = FALSE))

  expect_warning(conv_obj$solver <- "direct", "convolution")
  expect_equal(conv_obj$solver, "conv")
})
