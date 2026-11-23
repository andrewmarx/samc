# Tests for SSF models. An SSF model has the same state space as an RW model,
# so every metric should match the same metric calculated from a samc object
# built directly from the SSF model's P matrix.

ssf_res = terra::rast(baselist$res * masklist$mask2, extent = terra::ext(0, 8, 0, 9), crs = "local")
ssf_abs = terra::rast(baselist$abs * masklist$mask2, extent = terra::ext(0, 8, 0, 9), crs = "local")
ssf_fid = terra::rast(baselist$fid * masklist$mask2, extent = terra::ext(0, 8, 0, 9), crs = "local")
ssf_init = terra::rast(baselist$init * masklist$mask2, extent = terra::ext(0, 8, 0, 9), crs = "local")

ssf_origin = 5
ssf_dest = 40
ssf_time = 10

for (dir in c(4, 8)) {
  invisible(utils::capture.output({
    ssf_obj = samc(ssf_res, ssf_abs, ssf_fid,
                   model = list(name = "SSF", fun = function(x) 1 / x[2], dir = dir, sym = FALSE, ssc = 1))
    ref_obj = samc(ssf_obj$p_matrix)
  }))

  ssf_pv = terra::values(ssf_init, mat = FALSE)
  ssf_pv = ssf_pv[is.finite(ssf_pv)]

  # Calculate a metric with both objects and return the results as plain vectors
  ssf_cmp = function(fun, ..., init = FALSE) {
    if (init) {
      res = fun(ssf_obj, ..., init = ssf_init)
      ref = fun(ref_obj, ..., init = ssf_pv)
    } else {
      res = fun(ssf_obj, ...)
      # The P matrix object warns that its number of clumps is unknown
      ref = suppressWarnings(fun(ref_obj, ...))
    }

    expect_equal(as.vector(unlist(res)), as.vector(unlist(ref)), tolerance = 1e-8)
  }

  test_that(paste("SSF model does not over-allocate the Q matrix for dir =", dir), {
    q = ssf_obj@data@f

    expect_equal(length(q@i), q@p[length(q@p)])
    expect_equal(length(q@x), q@p[length(q@p)])
  })

  test_that(paste("SSF location inputs are validated like RW for dir =", dir), {
    expect_error(visitation(ssf_obj, origin = 0))
    expect_error(visitation(ssf_obj, origin = length(ssf_pv) + 1))
  })

  test_that(paste("SSF visitation() matches P matrix results for dir =", dir), {
    ssf_cmp(visitation, origin = ssf_origin)
    ssf_cmp(visitation, dest = ssf_dest)
    ssf_cmp(visitation, origin = ssf_origin, dest = ssf_dest)
    ssf_cmp(visitation, init = TRUE)
    ssf_cmp(visitation, origin = ssf_origin, time = ssf_time)
    ssf_cmp(visitation, dest = ssf_dest, time = ssf_time)
    ssf_cmp(visitation, time = ssf_time, init = TRUE)
  })

  test_that(paste("SSF dispersal() matches P matrix results for dir =", dir), {
    ssf_cmp(dispersal, dest = ssf_dest)
    ssf_cmp(dispersal, origin = ssf_origin, dest = ssf_dest)
    ssf_cmp(dispersal, dest = ssf_dest, time = ssf_time)
    ssf_cmp(dispersal, origin = ssf_origin, dest = ssf_dest, time = ssf_time)
    ssf_cmp(dispersal, dest = ssf_dest, init = TRUE)
    ssf_cmp(dispersal, dest = ssf_dest, time = ssf_time, init = TRUE)
  })

  test_that(paste("SSF distribution() matches P matrix results for dir =", dir), {
    ssf_cmp(distribution, origin = ssf_origin, time = ssf_time)
    ssf_cmp(distribution, dest = ssf_dest, time = ssf_time)
    ssf_cmp(distribution, origin = ssf_origin, dest = ssf_dest, time = ssf_time)
    ssf_cmp(distribution, time = ssf_time, init = TRUE)
  })

  test_that(paste("SSF mortality() matches P matrix results for dir =", dir), {
    ssf_cmp(mortality, origin = ssf_origin)
    ssf_cmp(mortality, dest = ssf_dest)
    ssf_cmp(mortality, origin = ssf_origin, dest = ssf_dest)
    ssf_cmp(mortality, init = TRUE)
    ssf_cmp(mortality, origin = ssf_origin, time = ssf_time)
    ssf_cmp(mortality, dest = ssf_dest, time = ssf_time)
    ssf_cmp(mortality, time = ssf_time, init = TRUE)
  })

  test_that(paste("SSF survival() matches P matrix results for dir =", dir), {
    ssf_cmp(survival)
    ssf_cmp(survival, init = TRUE)
  })

  test_that(paste("SSF cond_passage() matches P matrix results for dir =", dir), {
    ssf_cmp(cond_passage, dest = ssf_dest)
    ssf_cmp(cond_passage, origin = ssf_origin, dest = ssf_dest)
  })
}
