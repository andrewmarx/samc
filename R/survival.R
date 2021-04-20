# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate survival metrics
#'
#' Calculates the expected amount of time that individuals survive
#' in the landscape.
#'
#' \eqn{z=(I-Q)^{-1}{\cdot}1=F{\cdot}1}
#' \itemize{
#'   \item \strong{survival(samc)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. The value of element \emph{i} is the expected
#' amount of time that individuals survive when starting at location \emph{i}.
#' }
#'
#' \eqn{\psi^Tz}
#' \itemize{
#'   \item \strong{survival(samc, occ)}
#'
#' The result is a numeric that represents the expected time that any individual
#' stays in the landscape before death, regardless of the initial location.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-occ
#'
#' @return A \code{vector} or a \code{numeric}
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "survival",
  function(samc, occ) {
    standardGeneric("survival")
  })

# survival(samc) ----
#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", occ = "missing"),
  function(samc) {
    q <- samc$q_matrix
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    z = .f1(q)

    return(as.vector(z))
  })

# survival(samc, occ) ----
#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", occ = "RasterLayer"),
  function(samc, occ) {
    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    sv <- survival(samc)

    surv <- pv %*% sv

    return(as.numeric(surv))
  })

#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", occ = "matrix"),
  function(samc, occ) {
    occ <- .rasterize(occ)

    return(survival(samc, occ))
  })
