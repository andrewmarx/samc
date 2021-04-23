# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R
NULL


#' Calculate visitation metrics
#'
#' Calculates the number of times that transient states are visited before absorption.
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

# visitation(samc) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "missing", dest = "missing"),
  function(samc){
    if (!samc@override)
      stop("This version of the visitation() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    q <- samc$q_matrix
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1
    n <- Matrix::solve(q)
    return(as.matrix(n))
  })

# visitation(samc, origin) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "location", dest = "missing"),
  function(samc, origin){
    if (length(origin) != 1)
      stop("origin can only contain a single value for this version of the function", call. = FALSE)

    origin = .process_locations(samc, origin)

    q <- samc$q_matrix
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    r <- .f_row(q, origin);
    return(as.vector(r))
  })

# visitation(samc, dest) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "missing", dest = "location"),
  function(samc, dest){
    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest <- .process_locations(samc, dest)

    q <- samc$q_matrix
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    r <- .f_col(q, dest);
    return(as.vector(r))
  })

# visitation(samc, origin, dest) ----
#' @rdname visitation
setMethod(
  "visitation",
  signature(samc = "samc", origin = "location", dest = "location"),
  function(samc, origin, dest){
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
