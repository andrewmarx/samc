context("Convolution")


conv_samc = samc(testlist[[1]]$res,
                 pmax(testlist[[1]]$abs, 0.00001),
                 testlist[[1]]$fid,
                 model = list(fun = "1/mean(x)", dir = 8, sym = TRUE),
                 options = list(threads = 1,
                                method = "conv",
                                override = FALSE))

samc_obj = samc(testlist[[1]]$res,
                pmax(testlist[[1]]$abs, 0.00001),
                testlist[[1]]$fid,
                model = list(fun = function(x) 1/mean(x), dir = 8, sym = TRUE))


test_that("Convolution short-term distribution()", {
  r1 = distribution(samc_obj, testlist[[1]]$init, time = time)
  r2 = distribution(conv_samc, testlist[[1]]$init, time = time)

  expect_equal(r1, r2)
})

test_that("Convolution short-term mortality()", {
  r1 = mortality(samc_obj, testlist[[1]]$init, time = time)
  r2 = mortality(conv_samc, testlist[[1]]$init, time = time)

  expect_equal(r1, r2)
})

test_that("Convolution short-term visitation()", {
  r1 = as.vector(visitation(samc_obj, testlist[[1]]$init, time = time))
  r2 = visitation(conv_samc, testlist[[1]]$init, time = time)

  expect_equal(r1, r2)
})

test_that("Convolution long-term mortality()", {
  r1 = mortality(samc_obj, testlist[[1]]$init)
  r2 = mortality(conv_samc, testlist[[1]]$init)

  expect_equal(r1, r2)
})

test_that("Convolution long-term visitation()", {
  r1 = as.vector(visitation(samc_obj, testlist[[1]]$init))
  r2 = visitation(conv_samc, testlist[[1]]$init)

  expect_equal(r1, r2)
})
