# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate distribution metrics
#'
#' Calculate the probability of being at a transient state at a specific time.
#'
#' \eqn{Q^t}
#' \itemize{
#'   \item \strong{distribution(samc, time)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the probability of being
#' at transient state \eqn{\mathit{j}} after \eqn{\mathit{t}} time steps if starting
#' at transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{distribution(samc, origin, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability
#' of being at transient state \eqn{\mathit{j}} after \eqn{\mathit{t}} time steps
#' if starting at transient state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{distribution(samc, dest, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the probability
#' of being at transient state \eqn{\mathit{j}} after \eqn{\mathit{t}} time steps
#' if starting at transient state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{distribution(samc, origin, dest, time)}
#'
#' The result is a numeric value that is the probability of being at a transient
#' state \eqn{\mathit{j}} after \eqn{\mathit{t}} time steps if starting at transient
#' state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#' }
#'
#' \eqn{\psi^TQ^t}
#' \itemize{
#'   \item \strong{distribution(samc, init, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability
#' of being at transient state \eqn{\mathit{i}} after \eqn{\mathit{t}} time steps
#' given an initial state \eqn{\psi}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-init
#' @template param-origin
#' @template param-dest
#' @template param-time
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "distribution",
  function(samc, init, origin, dest, time) {
    standardGeneric("distribution")
  })

# distribution(samc, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the distribution() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer", call. = FALSE)

    q <- as.matrix(samc$q_matrix)

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
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "numeric"),
  function(samc, origin, time) {
    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1) stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    init = .process_locations(samc, origin)

    return(distribution(samc, init, time = time))
  })

# distribution(samc, dest, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "numeric"),
  function(samc, dest, time) {
    .disable_conv(samc)
    .disable_crw(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest <- .process_locations(samc, dest, map = FALSE)
    .validate_time_steps(time)

    q <- samc$q_matrix

    time <- c(1, time)

    mov <- .qpow_col(q, dest, time)
    mov <- lapply(mov, as.vector)

    if (length(mov) == 1) {
      return((mov[[1]]))
    } else {
      return(mov)
    }
  })

# distribution(samc, origin, dest, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "numeric"),
  function(samc, origin, dest, time) {
    .disable_conv(samc)
    .disable_crw(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest <- .process_locations(samc, dest, map = FALSE)

    mov <- distribution(samc, origin = origin, time = time)

    if (is.list(mov)){
      return(lapply(mov, "[", dest))
    } else if (is.vector(mov)) {
      return(mov[dest])
    } else {
      stop("This should not have been possible. Please submit a report with a fully reproducible and simplified example.", call. = FALSE)
    }
  })

# distribution(samc, init, time) ----
#' @rdname distribution
setMethod(
  "distribution",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, init, time) {
    check(samc, init)

    pv <- .process_init(samc, init)

    .validate_time_steps(time)

    if (samc@solver %in% c("direct", "iter")) {
      q = samc$q_matrix

      time = c(0, time)

      res = .qpow_row(q, pv, time)

      res = lapply(res, as.vector)

      if (length(res) == 1) {
        return(res[[1]])
      } else {
        return(res)
      }
    } else if (samc@solver == "conv") {
      results_list <- samc:::.convolution_short(time, samc@conv_cache, pv, samc@threads)

      res = as.vector(results_list$dist[[1]])

      return(res)
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })
