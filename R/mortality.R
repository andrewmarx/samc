# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R visitation.R
NULL


#' Calculate mortality metrics
#'
#' Calculates the probability of absorption at individual transient states.
#'
#' \eqn{\tilde{B}_t = (\sum_{n=0}^{t-1} Q^n) \tilde{R}}
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
#'   \item \strong{mortality(samc, occ, time)}
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
#'   \item \strong{mortality(samc, occ)}
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
#' @template param-occ
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
  function(samc, occ, origin, dest, time) {
    standardGeneric("mortality")
  })

# mortality(samc, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {
    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer", call. = FALSE)

    # TODO: remove as.matrix call, which is needed to convert from a sparse to
    # dense matrix for the %^% operator, which means removing expm as a dependency
    q <- as.matrix(samc$q_matrix)
    r <- matrix(0, nrow = nrow(q), ncol = nrow(q))
    diag(r) <- rowSums(samc$r_matrix)
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
  signature(samc = "samc", occ = "missing", origin = "location", dest = "missing", time = "numeric"),
  function(samc, origin, time) {
    if (length(origin) != 1)
      stop("origin can only contain a single location for this version of the function", call. = FALSE)

    origin <- .process_locations(samc, origin)
    .validate_time_steps(time)

    q <- samc$q_matrix

    rdg <- rowSums(samc$r_matrix)

    time <- c(1, time)

    mort <- .sum_qpow_row(q, origin, time)

    mort <- lapply(mort, function(x){as.vector(x * rdg)})

    if (length(mort) == 1) {
      return(mort[[1]])
    } else {
      return(mort)
    }
  })

# mortality(samc, dest, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "location", time = "numeric"),
  function(samc, dest, time) {
    if (length(dest) != 1)
      stop("dest can only contain a single location for this version of the function", call. = FALSE)

    dest <- .process_locations(samc, dest)
    .validate_time_steps(time)

    q <- samc$q_matrix

    rdg <- rowSums(samc$r_matrix)

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
  signature(samc = "samc", occ = "missing", origin = "location", dest = "location", time = "numeric"),
  function(samc, origin, dest, time) {
    dest <- .process_locations(samc, dest)

    mort <- mortality(samc, origin = origin, time = time)

    if (is.list(mort)){
      return(lapply(mort, "[", dest))
    } else if (is.vector(mort)) {
      return(mort[dest])
    } else {
      stop("This should not have been possible. Please submit a report with a fully reproducible and simplified example.", call. = FALSE)
    }
  })

# mortality(samc, occ, time) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, occ, time) {
    .validate_time_steps(time)

    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    q <- samc$q_matrix
    Rdiag <- rowSums(samc$r_matrix)

    time <- c(1, time)
    mort <- .sum_psiqpow(q, pv, time)

    mort <- lapply(mort, function(x){as.vector(x * Rdiag)})

    if (length(mort) == 1) {
      return(mort[[1]])
    } else {
      return(mort)
    }
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, occ, time) {
    occ <- .rasterize(occ)

    return(mortality(samc, occ, time = time))
  })

# mortality(samc) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc) {
    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nSee the documentation for details.", call. = FALSE)

    f <- visitation(samc)
    gc()
    dimnames(f) <- dimnames(samc$q_matrix) # Not sure why dimnames aren't carrying through later calculations

    rdg <- rowSums(samc$r_matrix)
    r <- Matrix::sparseMatrix(i = 1:length(rdg),
                              j = 1:length(rdg),
                              x = rdg,
                              index1 = TRUE)

    # TODO f %*% r can be simplified to an elementwise multiplication of the matrix columns by the corresponding elements in the rdg vector. This might be helpful for memory allocations and performance.
    mort <- f %*% r
    dimnames(mort) <- dimnames(samc$q_matrix) # See above dimnames comment

    gc()

    if (ncol(samc$r_matrix) > 1) {
      mort_list <- list()
      for (n in colnames(samc$r_matrix)) {
        Matrix::diag(r) <- samc$r_matrix[, n]
        mort_list[[n]] <- f %*% r
        dimnames(mort_list[[n]]) <- dimnames(samc$q_matrix) # See above dimnames comment
      }
      mort_list$total <- mort
      return(mort_list)
    } else {
      return(mort)
    }
  })

# mortality(samc, origin) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "location", dest = "missing", time = "missing"),
  function(samc, origin) {
    vis <- visitation(samc, origin = origin)
    names(vis) <- rownames(samc$q_matrix)

    rdg <- rowSums(samc$r_matrix)
    mort <- vis * rdg

    if (ncol(samc$r_matrix) > 1) {
      mort_list <- list()
      for (n in colnames(samc$r_matrix)) {
        mort_list[[n]] <- vis * samc$r_matrix[, n]
      }
      mort_list$total <- mort
      return(mort_list)
    } else {
      return(mort)
    }
  })

# mortality(samc, dest) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "location", time = "missing"),
  function(samc, dest) {
    vis <- visitation(samc, dest = dest)
    names(vis) <- rownames(samc$q_matrix)

    rdg <- rowSums(samc$r_matrix)
    mort <- vis * rdg[dest]

    if (ncol(samc$r_matrix) > 1) {
      mort_list <- list()
      for (n in colnames(samc$r_matrix)) {
        mort_list[[n]] <- vis * samc$r_matrix[dest, n]
      }
      mort_list$total <- mort
      return(mort_list)
    } else {
      return(mort)
    }
  })

# mortality(samc, origin, dest) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "location", dest = "location", time = "missing"),
  function(samc, origin, dest) {
    if(length(origin) != length(dest))
      stop("The 'origin' and 'dest' parameters must have the same number of values", call. = FALSE)

    origin <- .process_locations(samc, origin)
    dest <- .process_locations(samc, dest)

    rdg <- rowSums(samc$r_matrix)

    results <- vector(mode = "numeric", length = length(origin))

    for (d in unique(dest)) {
      vis <- visitation(samc, dest = d)
      results[dest == d] <- vis[origin[dest == d]]
    }
    names(results) <- rownames(samc$q_matrix)[dest]

    mort <- results * rdg[dest]

    if (ncol(samc$r_matrix) > 1) {
      mort_list <- list()
      for (n in colnames(samc$r_matrix)) {
        mort_list[[n]] <- results * samc$r_matrix[dest, n]
      }
      mort_list$total <- mort
      return(mort_list)
    } else {
      return(mort)
    }
  })

# mortality(samc, occ) ----
#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {
    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    q <- samc$q_matrix
    rdg <- rowSums(samc$r_matrix)

    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    pf <- .psif(q, pv)
    names(pf) <- rownames(samc$q_matrix)

    mort <- pf * rdg

    if (ncol(samc$r_matrix) > 1) {
      mort_list <- list()
      for (n in colnames(samc$r_matrix)) {
        mort_list[[n]] <- pf * samc$r_matrix[, n]
      }
      mort_list$total <- mort
      return(mort_list)
    } else {
      return(mort)
    }
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {
    occ <- .rasterize(occ)

    return(mortality(samc, occ))
  })
