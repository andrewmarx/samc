# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R visitation.R
NULL


#' Calculate dispersal metrics
#'
#' Calculates the probability of individuals visiting locations
#'
#'
#' \eqn{\tilde{D}_{jt}=(\sum_{n=0}^{t-1}\tilde{Q}^n)\tilde{q}_j}
#' \itemize{
#'   \item \strong{dispersal(samc, dest, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the probability
#' of visiting transient state \eqn{\mathit{j}} within \eqn{\mathit{t}} or fewer
#' time steps if starting at transient state \eqn{\mathit{i}}.
#'
#' Note: Given the current derivation, when \eqn{\mathit{i=j}}, then \eqn{\mathbf{v}_i}
#' is unknown and has been set to \code{NA}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#' }
#'
#' \eqn{\psi^T\tilde{D}_{jt}}
#' \itemize{
#'   \item \strong{dispersal(samc, init, dest, time)}
#'
#' The result is a numeric that is the probability of visiting transient state \eqn{\mathit{j}}
#' within \eqn{\mathit{t}} or fewer time steps given an initial state \eqn{\psi}
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#' }
#'
#' \eqn{D=(F-I)diag(F)^{-1}}
#' \itemize{
#'   \item \strong{dispersal(samc)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the probability of visiting
#' transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{dispersal(samc, origin)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability
#' of visiting transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{dispersal(samc, dest)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the probability
#' of visiting transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{dispersal(samc, origin, dest)}
#'
#' The result is a numeric value that is the probability of visiting transient
#' state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#' }
#'
#' \eqn{\psi^TD}
#' \itemize{
#'   \item \strong{dispersal(samc, init)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability
#' of visiting transient state \eqn{\mathit{j}} given an initial state \eqn{\psi}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{dispersal(samc, init, dest)}
#'
#' The result is a numeric value that is the probability of visiting transient
#' state \eqn{\mathit{j}} given an initial state \eqn{\psi}.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-init
#' @template param-origin
#' @template param-dest
#' @template param-time
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "dispersal",
  function(samc, init, origin, dest, time) {
    standardGeneric("dispersal")
  })

# dispersal(samc, origin, dest, time) ----
#' @rdname dispersal
setMethod(
  "dispersal", # TODO add unit tests
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "numeric"),
  function(samc, origin, dest, time) {
    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1)
        stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    origin = .process_locations(samc, origin)
    init = .map_location(samc, origin)

    return(dispersal(samc, init, dest=dest, time = time))
  })

# dispersal(samc, dest, time) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "numeric"),
  function(samc, dest, time) {
    .disable_conv(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest <- .process_locations(samc, dest)
    .validate_time_steps(time)

    if (samc@model$name == "RW") {
      vec = logical(samc@nodes)
      vec[dest] = TRUE
    } else if (samc@model$name == "CRW") {
      vec = (samc@crw_map[,1] == dest)
    } else {
      stop("Unexpected model", call. = FALSE)
    }

    q = samc$q_matrix

    qv = q[, vec, drop = FALSE]
    qv[vec ,] = 0
    qv = Matrix::rowSums(qv)

    q[, vec] = 0
    q[vec, ] = 0

    q2 = q
    q2@x = -q2@x
    Matrix::diag(q2) = Matrix::diag(q2) + 1

    time <- c(0, time)

    if (samc@solver == "iter") {
      res = .sum_qn_q_iter(q, q2, qv, time)
    } else {
      res = .sum_qn_q(q, q2, qv, time)
    }

    res = lapply(res, as.vector)

    if (samc@model$name == "CRW") {
      pv = samc@prob_mat
      pv = pv[!is.na(pv)]

      res = lapply(res, function(x) .summarize_crw(samc, pv * x, sum))
    }

    if (length(res) == 1) {
      return(res[[1]])
    } else {
      return(res)
    }
  })

# dispersal(samc, init, dest, time) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "location", time = "numeric"),
  function(samc, init, dest, time) {
    .disable_conv(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    check(samc, init)

    dest <- .process_locations(samc, dest)

    pv <- .process_init(samc, init)

    if (samc@model$name == "CRW") pv = .summarize_crw(samc, pv, sum)

    d <- dispersal(samc, dest = dest, time = time)

    pv <- pv[-dest]

    if (is.list(d)) {
      return(lapply(d, function(x){as.numeric(pv %*% x[-dest])}))
    } else {
      return(as.numeric(pv %*% d[-dest]))
    }
  })


# dispersal(samc) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc) {
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the dispersal() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    f <- visitation(samc)
    gc()
    fdg <- 1 / Matrix::diag(f)
    fdg_mat <- Matrix::sparseMatrix(i = 1:length(fdg),
                                    j = 1:length(fdg),
                                    x = fdg,
                                    index1 = TRUE)
    # TODO Check if 'diag(f) <- ' resulting in extra memory allocation.
    Matrix::diag(f) <- Matrix::diag(f) - 1
    gc()
    d_mat <- f %*% fdg_mat

    return(d_mat)
  })

# dispersal(samc, origin) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "missing"),
  function(samc, origin) {
    .disable_conv(samc)

    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1)
        stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    origin = .process_locations(samc, origin)
    init = .map_location(samc, origin)

    return(dispersal(samc, init))
  })

# dispersal(samc, dest) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "missing"),
  function(samc, dest) {
    .disable_conv(samc)

    dest <- .process_locations(samc, dest)

    f_col <- visitation(samc, dest = dest)
    fjj <- f_col[dest]
    f_col[dest] <- f_col[dest] - 1

    result <- as.vector(f_col/fjj)
    names(result) <- samc$names

    return(result)
  })

# dispersal(samc, origin, dest) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "missing"),
  function(samc, origin, dest) {
    .disable_conv(samc)

    origin <- .process_locations(samc, origin)
    dest <- .process_locations(samc, dest)

    if(length(origin) != length(dest))
      stop("The 'origin' and 'dest' parameters must have the same number of values", call. = FALSE)

    result <- vector(mode = "numeric", length = length(origin))

    for (d in unique(dest)) {
      # Using dispersal(samc, dest) because dispersal(samc, origin) is not optimized
      t <- dispersal(samc, dest = d)
      result[dest == d] <- t[origin[dest == d]]
    }

    return(result)
  })

# dispersal(samc, init) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "missing"),
  function(samc, init) {
    .disable_conv(samc)

    check(samc, init)

    if (!samc@.cache$dgf_exists) {
      if (samc@solver == "iter") {
        dg <- samc:::.diagf_par_iter(samc@data@f, samc@threads)
      } else {
        dg <- samc:::.diagf_par(samc@data@f, samc@threads)
      }

      samc@.cache$dgf <- dg
      samc@.cache$dgf_exists <- TRUE
    }

    vis = visitation(samc, init)

    dg = samc@.cache$dgf
    init = .process_init(samc, init)

    if (samc@model$name == "CRW") {
      dg = .summarize_crw(samc, dg, sum) - .summarize_crw(samc, dg, length) + 1
      init = .summarize_crw(samc, init, sum)
    }

    return((vis - init)/dg)
  })


# dispersal(samc, init, dest) ----
#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "location", time = "missing"),
  function(samc, init, dest) {
    .disable_conv(samc)

    check(samc, init)

    pv <- .process_init(samc, init)

    dj <- dispersal(samc, dest = dest)

    if (samc@model$name == "CRW") {
      pv = .summarize_crw(samc, pv, sum)
    }

    return(as.numeric(pv %*% dj))
  })
