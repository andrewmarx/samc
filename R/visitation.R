# Copyright (c) 2019-2023 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R
NULL


#' Calculate visitation metrics
#'
#' Calculates the number of times that transient states are visited before absorption.
#'
#'
#' \eqn{\tilde{F}_{t} = (\sum_{n=0}^{t-1}{Q}^n)}
#' \itemize{
#'   \item \strong{visitation(samc, time)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the number of times that
#' transient state \eqn{\mathit{j}} is visited after \eqn{\mathit{t}} time steps
#' if starting at transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{visitation(samc, origin, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the number
#' of times that transient state \eqn{\mathit{j}} is visited after \eqn{\mathit{t}}
#' time steps if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, dest, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the number
#' of times that transient state \eqn{\mathit{j}} is visited after \eqn{\mathit{t}}
#' time steps if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, origin, dest, time)}
#'
#' The result is a numeric value that is the number of times transient state
#' \eqn{\mathit{j}} is visited after \eqn{\mathit{t}} time steps if starting at
#' transient state \eqn{\mathit{i}}.
#' }
#'
#' \eqn{\psi^T \tilde{F}_{t}}
#' \itemize{
#'   \item \strong{visitation(samc, init, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the number
#' of times that transient state \eqn{\mathit{j}} is visited after \eqn{\mathit{t}}
#' time steps before absorption given an initial state \eqn{\psi}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, init, dest, time)}
#'
#' The result is a numeric value that is the number of times transient state
#' \eqn{\mathit{j}} is visited after \eqn{\mathit{t}} time steps given an initial
#' state \eqn{\psi}.
#' }
#'
#'
#' \eqn{F = (I-Q)^{-1}}
#' \itemize{
#'   \item \strong{visitation(samc)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the number of times that
#' transient state \eqn{\mathit{j}} is visited before absorption if starting at
#' transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{visitation(samc, origin)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the number
#' of times that transient state \eqn{\mathit{j}} is visited before absorption if
#' starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, dest)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the number
#' of times that transient state \eqn{\mathit{j}} is visited before absorption if
#' starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, origin, dest)}
#'
#' The result is a numeric value that is the number of times transient state
#' \eqn{\mathit{j}} is visited before absorption if starting at transient
#' state \eqn{\mathit{i}}.
#' }
#'
#' \eqn{\psi^TF}
#' \itemize{
#'   \item \strong{visitation(samc, init)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the number
#' of times that transient state \eqn{\mathit{j}} is visited before absorption
#' given an initial state \eqn{\psi}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{visitation(samc, init, dest)}
#'
#' The result is a numeric value that is the number of times transient state
#' \eqn{\mathit{j}} is visited before absorption given an initial state \eqn{\psi}.
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
  "visitation",
  function(samc, init, origin, dest, time) {
    standardGeneric("visitation")
  })

# TODO tests for short-term visitation metrics
# visitation(samc, time) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time){
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the visitation() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer", call. = FALSE)

    q <- as.matrix(samc$q_matrix)
    q2 = base::diag(nrow(q))
    res <- q2

    for (i in 1:(time-1)) {
      q2 = q2 %*% q
      res <- res + q2
    }

    return(res)
  })

# visitation(samc, origin, time) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "numeric"),
  function(samc, origin, time){

    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1) stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }
    .validate_time_steps(time)

    origin = .process_locations(samc, origin)
    init = .map_location(samc, origin)

    return(visitation(samc, init, time = time))
  })

# visitation(samc, dest, time) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "numeric"),
  function(samc, dest, time){
    .disable_conv(samc)
    .disable_crw(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest = .process_locations(samc, dest)
    .validate_time_steps(time)

    q = samc$q_matrix

    time = c(1, time)

    ft = .sum_qpow_col(q, dest, time)

    if (length(ft) == 1) {
      return(ft[[1]])
    } else {
      return(ft)
    }
  })

# visitation(samc, origin, dest, time) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "numeric"),
  function(samc, origin, dest, time){
    .disable_conv(samc)
    .disable_crw(samc)

    dest = .process_locations(samc, dest)

    ft <- visitation(samc, origin = origin, time = time)

    if (is.list(ft)){
      return(lapply(ft, "[", dest))
    } else if (is.vector(ft)) {
      return(ft[dest])
    } else {
      stop("This should not have been possible. Please submit a report with a fully reproducible and simplified example.", call. = FALSE)
    }
  })

# visitation(samc, init, time) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, init, time){
    .validate_time_steps(time)

    check(samc, init)

    pv <- .process_init(samc, init)

    if (samc@solver %in% c("direct", "iter")) {
      q <- samc$q_matrix

      time <- c(1, time)
      ft <- .sum_qpow_row(q, pv, time)

      if (samc@model$name == "CRW") lapply(ft, function(x) .summarize_crw(samc, x, sum))

      if (length(ft) == 1) {
        return(ft[[1]])
      } else {
        return(ft)
      }
    } else if (samc@solver == "conv") {
      results_list = samc:::.convolution_short(time, samc@conv_cache, pv, samc@threads)

      res = as.vector(results_list$vis[[1]])

      return(res)
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })


# visitation(samc) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc){
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the visitation() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    n <- Matrix::solve(samc@data@f)
    return(as.matrix(n))
  })

# visitation(samc, origin) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "missing"),
  function(samc, origin){

    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1)
        stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    origin = .process_locations(samc, origin)
    init = .map_location(samc, origin)

    return(visitation(samc, init = init))
  })

# visitation(samc, dest) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "missing"),
  function(samc, dest){
    .disable_conv(samc)
    .disable_crw(samc)

    if (is.matrix(dest)) {

    } else {
      if (length(dest) != 1)
        stop("dest can only contain a single location for this version of the function", call. = FALSE)
    }

    dest <- .process_locations(samc, dest)

    if (samc@solver == "iter") {
      r <- .f_col_iter(samc@data@f, dest);
    } else {
      r <- .f_col(samc@data@f, dest, samc@.cache$sc);
    }

    return(as.vector(r))
  })

# visitation(samc, origin, dest) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "missing"),
  function(samc, origin, dest){
    .disable_conv(samc)
    .disable_crw(samc)

    origin <- .process_locations(samc, origin)
    dest <- .process_locations(samc, dest)

    if(length(origin) != length(dest))
      stop("The 'origin' and 'dest' parameters must have the same number of values", call. = FALSE)

    result <- vector(mode = "numeric", length = length(origin))

    for (o in unique(origin)) {
      # Using visitiation(samc, origin) because visitation(samc, dest) involves an extra transpose operation
      t <- visitation(samc, origin = o)
      result[origin == o] <- t[dest[origin == o]]
    }

    return(result)
  })

# visitation(samc, init) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "missing"),
  function(samc, init){
    check(samc, init)

    pv <- .process_init(samc, init)

    if (samc@solver %in% c("direct", "iter")) {
      if (samc@solver == "iter") {
        r <- .f_row_iter(samc@data@f, pv)
      } else {
        r <- .f_row(samc@data@f, pv, samc@.cache$sc)
      }

      r = as.vector(r)

      if (samc@model$name == "CRW") r = .summarize_crw(samc, r, sum)

      return(r)
    } else if (samc@solver == "conv") {

      results_list = samc:::.convolution_long(samc@conv_cache, pv, samc@threads)

      return(results_list$vis)
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })


# visitation(samc, init, dest) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "location", time = "missing"),
  function(samc, init, dest){
    .disable_crw(samc)
    .disable_conv(samc)

    check(samc, init)

    pv <- .process_init(samc, init)

    fj <- visitation(samc, dest = dest)

    return(as.numeric(pv %*% fj))
  })


#' Calculate net visitation
#'
#' Calculates the net number of times that transient states are visited before absorption.
#'
#' The \code{\link{visitation}} function calculates the
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-init
#' @template param-origin
#' @template param-dest
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "visitation_net",
  function(samc, init, origin, dest) {
    standardGeneric("visitation_net")
  })

# visitation_net(samc, origin, dest) ----
#' @rdname visitation_net
setMethod(
  "visitation_net",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location"),
  function(samc, origin, dest) {
    if (length(origin) != 1)
      stop("origin can only contain a single location for this version of the function", call. = FALSE)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    vis = visitation(samc, origin = origin)

    vq = vis * samc@data@f

    n_net = abs(Matrix::skewpart(vq))
    visit_net = as.vector(Matrix::colSums(n_net))
    visit_net[c(origin, dest)] = 2 * visit_net[c(origin, dest)]

    return(visit_net)
  })
