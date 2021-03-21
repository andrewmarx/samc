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

# distribution(samc, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {

    if (!samc@override)
      stop("This version of the distribution() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.", call. = FALSE)

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer", call. = FALSE)

    q <- as.matrix(samc@p[-nrow(samc@p), -nrow(samc@p)])

    res <- base::diag(nrow(q))

    for (i in 1:time) {
      res <- res %*% q
    }

    return(res)
  })

# distribution(samc, origin, time) ----
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
  signature(samc = "samc", occ = "missing", origin = "character", dest = "missing", time = "numeric"),
  function(samc, origin, time) {
    if (length(origin) != 1)
      stop("origin can only contain a single location for this version of the function", call. = FALSE)

    row_names <- rownames(samc@p)
    .validate_names(row_names[-length(row_names)], origin)

    return(distribution(samc, origin = match(origin, row_names), time = time))
  })

# distribution(samc, dest, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, dest, time) {
    .validate_locations(samc, dest)
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
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "character", time = "numeric"),
  function(samc, dest, time) {
    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    col_names <- colnames(samc@p)
    .validate_names(col_names[-length(col_names)], dest)

    return(distribution(samc, dest = match(dest, col_names), time = time))
  })

# distribution(samc, origin, dest, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "numeric", time = "numeric"),
  function(samc, origin, dest, time) {
    .validate_locations(samc, dest)

    mov <- distribution(samc, origin = origin, time = time)

    if (is.list(mov)){
      return(lapply(mov, "[", dest))
    } else if (is.vector(mov)) {
      return(mov[dest])
    } else {
      stop("This should not have been possible. Please submit a report with a fully reproducible and simplified example.", call. = FALSE)
    }
  })

#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", occ = "missing", origin = "character", dest = "character", time = "numeric"),
  function(samc, origin, dest, time) {
    row_names <- rownames(samc@p)
    .validate_names(row_names[-length(row_names)], origin)

    col_names <- colnames(samc@p)
    .validate_names(col_names[-length(col_names)], dest)

    return(distribution(samc,
                        origin = match(origin, row_names),
                        dest = match(dest, col_names),
                        time = time))
  })


#
# \psiQ^t
#

# distribution(samc, occ, time) ----
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

    occ <- .rasterize(occ)

    return(distribution(samc, occ = occ, time = time))
  })
