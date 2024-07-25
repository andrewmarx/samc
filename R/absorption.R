# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R
NULL


#' Calculate absorption metrics
#'
#' Calculates the probability of absorption for absorbing states rather
#' than individual transient states. This is distint from, yet very closely linked
#' to, the mortality() metric, which calculates the probability of absorption at
#' individual transient states. If the results of the mortality() metric are decomposed
#' into individual results for each absorbing state, then the sums of the individual
#' results for every transient state are equivalent to the results of the absorption()
#' metric.
#'
#' \eqn{A = F R}
#' \itemize{
#'   \item \strong{absorption(samc)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} if starting
#' at transient state \eqn{\mathit{i}}.
#'
#'   \item \strong{absorption(samc, origin)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_{k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} if starting
#' at transient state \eqn{\mathit{i}}.
#' }
#'
#' \eqn{\psi^T A}
#' \itemize{
#'   \item \strong{absorption(samc, init)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_{k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} given an
#' initial state \eqn{\psi}.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-init
#' @template param-origin
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "absorption",
  function(samc, init, origin) {
    standardGeneric("absorption")
  })

# TODO Add unit tests

# absorption(samc) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", init = "missing", origin = "missing"),
  function(samc) {
    .disable_conv(samc)

    c_abs = samc@data@c_abs
    if (any(dim(c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    if (samc@model$name == "CRW") {
      c_abs = apply(c_abs, 2, function(x) { x[samc@crw_map[, 1]] })
    }

    # TODO: possibly optimize
    abs_mat <- Matrix::solve(samc@data@f, c_abs)

    abs_mat <- as.matrix(abs_mat)

    if (samc@model$name == "CRW") {
      vec = as.vector(samc@prob_mat)
      vec = vec[!is.na(vec)]

      abs_mat = vec * abs_mat

      abs_mat = apply(abs_mat, 2, function(x) { .summarize_crw(samc, x, sum) })
    }

    colnames(abs_mat) <- colnames(samc@data@c_abs)

    return(abs_mat)
  })

# absorption(samc, origin) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", init = "missing", origin = "location"),
  function(samc, origin) {
    .disable_conv(samc)

    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1)
        stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    vis <- visitation(samc, origin = origin)

    result <- as.vector(vis %*% samc@data@c_abs)

    names(result) <- colnames(samc@data@c_abs)

    return(result)
  })

# absorption(samc, init) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", init = "ANY", origin = "missing"),
  function(samc, init) {
    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    pf = visitation(samc, init)

    result = as.vector(pf %*% samc@data@c_abs)
    names(result) = colnames(samc@data@c_abs)

    return(result)
  })
