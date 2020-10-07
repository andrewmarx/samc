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
#' amount of time that individuals survive when starting at location \emph{j}.
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

#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", occ = "missing"),
  function(samc) {
    q = samc@p[-nrow(samc@p),-nrow(samc@p)]
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    z = .f1(q)

    return(as.vector(z))
  })

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

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(survival(samc, occ))
  })
