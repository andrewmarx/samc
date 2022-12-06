# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

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
#'   \item \strong{absorption(samc, occ)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_{k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} given an
#' initial state \eqn{\psi}.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-occ
#' @template param-origin
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "absorption",
  function(samc, occ, origin) {
    standardGeneric("absorption")
  })


# absorption(samc) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "missing", origin = "missing"),
  function(samc) {
    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    # TODO: possibly optimize using C++
    abs_mat <- Matrix::solve(samc@data@f, samc@data@c_abs)

    abs_mat <- as.matrix(abs_mat)

    colnames(abs_mat) <- colnames(samc@data@c_abs)

    return(abs_mat)
  })

# absorption(samc, origin) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "missing", origin = "location"),
  function(samc, origin) {
    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    vis <- visitation(samc, origin = origin)

    result <- as.vector(vis %*% samc@data@c_abs)

    names(result) <- colnames(samc@data@c_abs)

    return(result)
  })

# absorption(samc, occ) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing"),
  function(samc, occ) {
    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    pf <- samc:::.psif(samc@data@f, pv)

    result <- as.vector(pf %*% samc@data@c_abs)
    names(result) <- colnames(samc@data@c_abs)

    return(result)
  })

#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "matrix", origin = "missing"),
  function(samc, occ) {
    if (any(dim(samc@data@c_abs) == 0)) stop("No absorption components defined in the samc object", call. = FALSE)

    occ <- samc:::.rasterize(occ)

    return(absorption(samc, occ))
  })
