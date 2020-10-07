# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate distribution metrics
#'
#' Calculate the probability of finding an individual at a given location at a
#' specific time.
#'
#' \eqn{Q^t}
#' \itemize{
#'   \item \strong{distribution(samc, time)}
#'
#' The result is a matrix where element (i,j) is the probability of being
#' at location j after t time steps if starting at location i.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use.
#'
#'   \item \strong{distribution(samc, origin, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple time steps)
#' where element j is the probability of being at location j after t time steps if
#' starting at a given origin.
#'
#'   \item \strong{distribution(samc, dest, time)}
#'
#' The result is a vector where element i is the probability of being
#' at a given destination after t time steps if starting at location i.
#'
#'   \item \strong{distribution(samc, origin, dest, time)}
#'
#' The result is a numeric value (single time step) or a list of numeric values
#' (multiple time steps) that is the probability of being at a given
#' destination after t time steps when beginning at a given origin.
#' }
#'
#' \eqn{\psi^TQ^t}
#' \itemize{
#'   \item \strong{distribution(samc, occ, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple
#' time steps) where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element \emph{i} is the unconditional
#' probability of finding an individual (or expected number of individuals) in
#' location \emph{i} after \emph{t} time steps.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-occ
#' @template param-origin
#' @template param-dest
#' @template param-time
#'
#' @return A \code{vector}
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "distribution",
  function(samc, occ, origin, dest, time) {
    standardGeneric("distribution")
  })

#
# Q^t
#

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {

    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.")

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer")

    q <- as.matrix(samc@p[-nrow(samc@p), -nrow(samc@p)])

    res <- base::diag(nrow(q))

    for (i in 1:time) {
      res <- res %*% q
    }

    return(res)
  })

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "missing", time = "numeric"),
  function(samc, origin, time) {

    .validate_time_steps(time)

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]

    time <- c(1, time)

    mov <- .qpow_row(q, origin, time)
    mov <- lapply(mov, as.vector)
    if (length(mov) == 1) {
      return((mov[[1]]))
    } else {
      return(mov)
    }
  })

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, dest, time) {

    .validate_time_steps(time)

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]

    time <- c(1, time)

    mov <- .qpow_col(q, dest, time)
    mov <- lapply(mov, as.vector)

    if (length(mov) == 1) {
      return((mov[[1]]))
    } else {
      return(mov)
    }
  })

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "numeric", time = "numeric"),
  function(samc, origin, dest, time) {

    mov <- distribution(samc, origin = origin, time = time)

    if (is.list(mov)){
      return(lapply(mov, "[", dest))
    } else if (is.vector(mov)) {
      return(mov[dest])
    } else {
      stop("Fatal error: This should not have been possible. Please submit a report with a fully reproducible and simplified example.")
    }
  })


#
# \psiQ^t
#

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, occ, time) {

    check(samc, occ)

    .validate_time_steps(time)

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    time <- c(0, time)

    res <- .psiq(q, pv, time)

    res <- lapply(res, as.vector)

    if (length(res) == 1) {
      return(res[[1]])
    } else {
      return(res)
    }
  })

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, occ, time) {

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(distribution(samc, occ = occ, time = time))
  })
