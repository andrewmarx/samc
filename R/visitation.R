# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL


#' Calculate visitation metrics
#'
#' Calculates the number of times that individuals from each location visit each
#' location in the landscape before death.
#'
#'
#' \eqn{F = (I-Q)^{-1}}
#' \itemize{
#'   \item \strong{visitation(samc)}
#'
#' The result is a matrix where element (\emph{i},\emph{j}) is the expected
#' number of times an individual that starts in \emph{i} uses \emph{j} before it
#' dies.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use.
#'
#'   \item \strong{visitation(samc, origin)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element \emph{j} is the number of times that an
#' individual starting at the origin visits location \emph{j} before it dies.
#'
#'   \item \strong{visitation(samc, dest)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element \emph{i} is the number of times that an
#' individual starting at location \emph{i} visits the destination before it
#' dies.
#'
#'   \item \strong{visitation(samc, origin, dest)}
#'
#' The result is a numeric value that is the expected number of times an
#' individual starting at the origin visits the destination before it dies.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-origin
#' @template param-dest
#'
#' @return A matrix, a vector, or a numeric
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "visitation",
  function(samc, origin, dest) {
    standardGeneric("visitation")
  })

#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "missing", dest = "missing"),
  function(samc){

    if (!samc@override)
      stop("This version of the visitation() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.")

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1
    n <- Matrix::solve(q)
    return(as.matrix(n))
  })

#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "numeric", dest = "missing"),
  function(samc, origin){
    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    r <- .f_row(q, origin);
    return(as.vector(r))
  })

#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "missing", dest = "numeric"),
  function(samc, dest){
    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    r <- .f_col(q, dest);
    return(as.vector(r))
  })

#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "numeric", dest = "numeric"),
  function(samc, origin, dest){

    v <- visitation(samc, origin)

    return(v[dest])
  })
