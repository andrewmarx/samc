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
#' \eqn{B = F R}
#' \itemize{
#'   \item \strong{absorption(samc)}
#'
#' The result is a matrix \eqn{\mathbf{M}} where \eqn{\mathbf{M}_{i,k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} if starting
#' at transient state \eqn{\mathit{i}}.
#'
#'   \item \strong{absorption(samc, origin)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_{k}} is the
#' probability of absorption due to absorbing state \eqn{\mathit{k}} if starting
#' at transient state \eqn{\mathit{i}}.
#'
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

    q <- samc$q_matrix

    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    # TODO: possibly optimize using C++
    abs_mat <- Matrix::solve(q, samc$r_matrix)

    abs_mat <- as.matrix(abs_mat)

    colnames(abs_mat) <- colnames(samc$r_matrix)

    return(abs_mat)
  })

# absorption(samc, origin) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "missing", origin = "location"),
  function(samc, origin) {
    vis <- visitation(samc, origin = origin)

    result <- as.vector(vis %*% samc$r_matrix)

    names(result) <- colnames(samc$r_matrix)

    return(result)
  })

# absorption(samc, occ) ----
#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing"),
  function(samc, occ) {
    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    q <- samc$q_matrix

    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    pf <- samc:::.psif(q, pv)

    result <- as.vector(pf %*% samc$r_matrix)
    names(result) <- colnames(samc$r_matrix)

    return(result)
  })

#' @rdname absorption
setMethod(
  "absorption",
  signature(samc = "samc", occ = "matrix", origin = "missing"),
  function(samc, occ) {
    occ <- samc:::.rasterize(occ)

    return(absorption(samc, occ))
  })
