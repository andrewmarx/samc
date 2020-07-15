# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL

#' Conditional Mean First Passage Time
#'
#' Calculate the mean number of steps to first passage
#'
#'
#' \eqn{\tilde{t}=\tilde{D}_j^{-1}\tilde{F}\tilde{D}_j{\cdot}1}
#' \itemize{
#'   \item \strong{cond_passage(samc, dest)}
#'
#' The result is a vector where each element corresponds to a cell in the landscape,
#' and can be mapped back to the landscape using the \code{\link{map}} function.
#' Element \emph{i} is the mean number of steps for the first passage time from
#' location \emph{i} conditional on absorption into \emph{j}
#'
#'   \item \strong{cond_passage(samc, origin, dest)}
#'
#' The result is a numeric value representing the mean number of steps for the
#' first passage time from a given origin conditional on absorption into a given
#' destination.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-origin
#' @template param-dest
#'
#' @return A numeric vector or a single numeric value
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "cond_passage",
  function(samc, origin, dest) {
    standardGeneric("cond_passage")
  })

#' @rdname cond_passage
setMethod(
  "cond_passage",
  signature(samc = "samc", origin = "missing", dest = "numeric"),
  function(samc, dest) {
    if (samc@source != "matrix") stop("This version of cond_passage() can only be used with samc-class objects created directly from a P matrix")

    if (dest %% 1 != 0 || dest < 1 || dest > (ncol(samc@p) - 1))
      stop("dest must be an integer that refers to a cell in the landscape")

    Q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    qj <- Q[-dest, dest]
    Qj <- Q[-dest, -dest]

    Qj@x <- -Qj@x
    Matrix::diag(Qj) <- Matrix::diag(Qj) + 1

    t <- .cond_t(Qj, qj)

    return(as.numeric(t))
  })

#' @rdname cond_passage
setMethod(
  "cond_passage",
  signature(samc = "samc", origin = "numeric", dest = "numeric"),
  function(samc, origin, dest) {

    if (origin %% 1 != 0 || origin < 1 || origin > (nrow(samc@p) - 1))
      stop("origin must be an integer that refers to a cell in the landscape")

    t <- cond_passage(samc, dest = dest)

    return(t[origin])
  })
