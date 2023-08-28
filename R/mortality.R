# Copyright (c) 2019-2023 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R visitation.R
NULL


#' Calculate mortality metrics
#'
#' Calculates the probability of absorption at individual transient states.
#'
#' \eqn{\tilde{B}_t = \tilde{F} \tilde{R}}
#' \itemize{
#'   \item \strong{mortality(samc, time)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the
#' probability of absorption at transient state \eqn{\mathit{j}} within \eqn{\mathit{t}}
#' or fewer steps if starting at transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{mortality(samc, origin, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability
#' of absorption at transient state \eqn{\mathit{j}} within \eqn{\mathit{t}} or
#' fewer steps if starting at transient state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{mortality(samc, dest, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the probability
#' of absorption at transient state \eqn{\mathit{j}} within \eqn{\mathit{t}} or
#' fewer steps if starting at transient state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{mortality(samc, origin, dest, time)}
#'
#' The result is a numeric value that is the probability of absorption at transient
#' state \eqn{\mathit{j}} within \eqn{\mathit{t}} or fewer time steps if starting
#' at transient state \eqn{\mathit{i}}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a numeric value for each time step.
#' }
#'
#' \eqn{\psi^T \tilde{B}_t}
#' \itemize{
#'   \item \strong{mortality(samc, init, time)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the unconditional
#' probability of absorption at transient state \eqn{\mathit{j}} within \eqn{\mathit{t}}
#' or fewer steps given an initial state \eqn{\psi}.
#'
#' If multiple time steps were provided as a vector, then the result will be an
#' ordered named list containing a vector for each time step.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#' }
#'
#' \eqn{B = F \tilde{R}}
#' \itemize{
#'   \item \strong{mortality(samc)}
#'
#' The result is a matrix \eqn{M} where \eqn{M_{i,j}} is the
#' probability of absorption at transient state \eqn{\mathit{j}} if starting at
#' transient state \eqn{\mathit{i}}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use (see \code{\link{samc-class}}).
#'
#'   \item \strong{mortality(samc, origin)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the probability of absorption
#' at transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{mortality(samc, dest)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the probability of absorption
#' at transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#'
#'   \item \strong{mortality(samc, origin, dest)}
#'
#' The result is a numeric value that is the probability of absorption
#' at transient state \eqn{\mathit{j}} if starting at transient state \eqn{\mathit{i}}.
#' }
#'
#' \eqn{\psi^T B}
#' \itemize{
#'   \item \strong{mortality(samc, init)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_j} is the unconditional
#' probability of absorption at transient state \eqn{\mathit{j}} given an initial
#' state \eqn{\psi}.
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
  "mortality",
  function(samc, init, origin, dest, time) {
    standardGeneric("mortality")
  })

# mortality(samc, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer", call. = FALSE)

    q <- as.matrix(samc$q_matrix)
    r <- matrix(0, nrow = nrow(q), ncol = nrow(q))
    diag(r) <- samc@data@t_abs
    qi <- diag(dim(q)[2])

    # Sum of geometric series
    # for t-1; so if t = 3, applies to t=2
    qt <- diag(nrow(q))

    for (i in 1:time) {
      qt <- qt %*% q
    }

    q_n <- solve(qi - q) %*% (qi - qt)

    bt <- q_n %*% r
    return(bt)
  })

# mortality(samc, origin, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "numeric"),
  function(samc, origin, time) {
    .disable_conv(samc)

    mort = visitation(samc, origin = origin, time = time)

    rdg <- samc@data@t_abs

    if (is.list(mort)) {
      return(lapply(mort, function(x){as.vector(x * rdg)}))
    } else {
      return(as.vector(mort * rdg))
    }
  })

# mortality(samc, dest, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "numeric"),
  function(samc, dest, time) {
    .disable_conv(samc)
    .disable_crw(samc)

    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest = .process_locations(samc, dest, map = FALSE)
    .validate_time_steps(time)

    q <- samc$q_matrix

    rdg <- samc@data@t_abs

    rdg[-dest] <- 0

    time <- c(1, time)

    mort <- .sum_qpowrv(q, rdg, time)

    mort <- lapply(mort, as.vector)

    if (length(mort) == 1) {
      return(mort[[1]])
    } else {
      return(mort)
    }
  })

# mortality(samc, origin, dest, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "numeric"),
  function(samc, origin, dest, time) {
    .disable_conv(samc)
    .disable_crw(samc)

    dest <- .process_locations(samc, dest, map = FALSE)

    mort <- mortality(samc, origin = origin, time = time)

    if (is.list(mort)){
      return(lapply(mort, "[", dest))
    } else if (is.vector(mort)) {
      return(mort[dest])
    } else {
      stop("This should not have been possible. Please submit a report with a fully reproducible and simplified example.", call. = FALSE)
    }
  })

# mortality(samc, init, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, init, time) {
    .disable_crw(samc)

    if (samc@solver %in% c("direct", "iter")) {
      mort = visitation(samc, init = init, time = time)

      rdg <- samc@data@t_abs

      if (is.list(mort)) {
        return(lapply(mort, function(x){as.vector(x * rdg)}))
      } else {
        return(as.vector(mort * rdg))
      }
    } else if (samc@solver == "conv") {

      res = visitation(samc, init, time = time)


      return(res * samc@data@t_abs)
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })

# mortality(samc) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc) {
    .disable_conv(samc)

    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    f <- visitation(samc)
    gc()
    dimnames(f) <- dimnames(samc$q_matrix) # Not sure why dimnames aren't carrying through later calculations

    rdg <- samc@data@t_abs
    r <- Matrix::sparseMatrix(i = 1:length(rdg),
                              j = 1:length(rdg),
                              x = rdg,
                              index1 = TRUE)

    # TODO f %*% r can be simplified to an elementwise multiplication of the matrix columns by the corresponding elements in the rdg vector. This might be helpful for memory allocations and performance.
    mort <- f %*% r
    dimnames(mort) <- dimnames(samc$q_matrix) # See above dimnames comment

    if (ncol(samc@data@c_abs) > 0) {
      mort <- list(total = mort)

      for (n in colnames(samc@data@c_abs)) {
        Matrix::diag(r) <- samc@data@c_abs[, n]
        # TODO f %*% r can be simplified to an elementwise multiplication of the matrix columns by the corresponding elements in the rdg vector. This might be helpful for memory allocations and performance.
        mort[[n]] <- f %*% r
        dimnames(mort[[n]]) <- dimnames(samc$q_matrix) # See above dimnames comment
        gc()
      }
    }

    return(mort)
  })

# mortality(samc, origin) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "location", dest = "missing", time = "missing"),
  function(samc, origin) {
    .disable_conv(samc)

    vis <- visitation(samc, origin = origin)
    names(vis) <- samc$names

    mort <- vis * samc@data@t_abs

    if (ncol(samc@data@c_abs) > 0) {
      mort <- list(total = mort)
      for (n in colnames(samc@data@c_abs)) {
        mort[[n]] <- vis * samc@data@c_abs[, n]
      }
    }

    return(mort)
  })

# mortality(samc, dest) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "missing", dest = "location", time = "missing"),
  function(samc, dest) {
    .disable_conv(samc)
    .disable_crw(samc)

    dest <- .process_locations(samc, dest, map = FALSE)

    vis <- visitation(samc, dest = dest)
    names(vis) <- samc$names

    mort <- vis * samc@data@t_abs[dest]

    if (ncol(samc@data@c_abs) > 0) {
      mort <- list(total = mort)
      for (n in colnames(samc@data@c_abs)) {
        mort[[n]] <- vis * samc@data@c_abs[dest, n]
      }
    }

    return(mort)
  })

# mortality(samc, origin, dest) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "missing", origin = "location", dest = "location", time = "missing"),
  function(samc, origin, dest) {
    .disable_conv(samc)
    .disable_crw(samc)

    if(length(origin) != length(dest))
      stop("The 'origin' and 'dest' parameters must have the same number of values", call. = FALSE)

    origin <- .process_locations(samc, origin, map = FALSE)
    dest <- .process_locations(samc, dest, map = FALSE)

    results <- vector(mode = "numeric", length = length(origin))

    for (d in unique(dest)) {
      vis <- visitation(samc, dest = d)
      results[dest == d] <- vis[origin[dest == d]]
    }
    names(results) <- samc$names[dest]


    mort <- results * samc@data@t_abs[dest]

    if (ncol(samc@data@c_abs) > 0) {
      mort <- list(total = mort)
      for (n in colnames(samc@data@c_abs)) {
        mort[[n]] <- results * samc@data@c_abs[dest, n]
      }
    }

    return(mort)
  })

# mortality(samc, init) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", init = "ANY", origin = "missing", dest = "missing", time = "missing"),
  function(samc, init) {
    .disable_crw(samc)

    check(samc, init)

    pv <- .process_init(samc, init)

    if (samc@solver %in% c("direct", "iter")) {
      if (samc@solver == "iter") {
        pf <- .f_row_iter(samc@data@f, pv)
      } else {
        pf <- .f_row(samc@data@f, pv, samc@.cache$sc)
      }

      names(pf) <- samc$names

      mort <- pf * samc@data@t_abs

      if (ncol(samc@data@c_abs) > 0) {
        mort <- list(total = mort)
        for (n in colnames(samc@data@c_abs)) {
          mort[[n]] <- pf * samc@data@c_abs[, n]
        }
      }

      return(mort)
    } else if (samc@solver == "conv") {

      res = visitation(samc, init)

      return(res *  samc@data@t_abs)
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })
