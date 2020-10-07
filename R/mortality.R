# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate mortality metrics
#'
#' Calculates the probability of experiencing mortality at specific locations.
#'
#' \eqn{\tilde{B}_t = (\sum_{n=0}^{t-1} Q^n) \tilde{R}}
#' \itemize{
#'   \item \strong{mortality(samc, time)}
#'
#' The result is a matrix where element (i,j) is the probability of experiencing
#' mortality at location j within t or fewer steps if starting at location i.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use.
#'
#'   \item \strong{mortality(samc, origin, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple
#' time steps) where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element j is the probability of experiencing
#' mortality at location j within t or fewer steps if starting at a given origin.
#'
#'   \item \strong{mortality(samc, dest, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple
#' time steps) where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element i is the probability of experiencing
#' mortality at a given destination within t or fewer steps if starting at
#' location i.
#'
#'   \item \strong{mortality(samc, origin, dest, time)}
#'
#' The result is a numeric value (single time step) or a list of numeric
#' values (multiple time steps) that is the probability of experiencing
#' mortality at a given destination within t or fewer steps if starting at a
#' given origin.
#' }
#'
#' \eqn{\psi^T \tilde{B}_t}
#' \itemize{
#'   \item \strong{mortality(samc, occ, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple
#' time steps) where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element j is the unconditional probability of
#' experiencing mortality at location j within t or fewer time steps.
#' }
#'
#' \eqn{B = F \tilde{R}}
#' \itemize{
#'   \item \strong{mortality(samc)}
#'
#' The result is a matrix where element (i,j) is the probability of experiencing
#' mortality at location j if starting at location i.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use.
#'
#'   \item \strong{mortality(samc, origin)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element j is the probability of experiencing
#' mortality at location j if starting at a given origin.
#'
#'   \item \strong{mortality(samc, dest)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element i is the probability of experiencing
#' mortality at a given destination if starting at location i.
#'
#'   \item \strong{mortality(samc, origin, dest)}
#'
#' The result is a numeric value that is the probability of experiencing
#' mortality at a given destination if starting at a given origin
#' }
#'
#' \eqn{\psi^T B}
#' \itemize{
#'   \item \strong{mortality(samc, occ)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element j is the unconditional probability of
#' experiencing mortality at location j, regardless of the initial state.
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
#' @return A matrix, vector, or numeric
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "mortality",
  function(samc, occ, origin, dest, time) {
    standardGeneric("mortality")
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, time) {

    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.")

    if (time %% 1 != 0 || time < 1 || length(time) > 1)
      stop("The time argument must be a single positive integer")

    # TODO: remove as.matrix call, which is needed to convert from a sparse to
    # dense matrix for the %^% operator, which means removing expm as a dependency
    q <- as.matrix(samc@p[-nrow(samc@p), -nrow(samc@p)])
    r <- matrix(0, nrow = nrow(q), ncol = nrow(q))
    diag(r) <- samc@p[-nrow(samc@p), ncol(samc@p)]
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

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "missing", time = "numeric"),
  function(samc, origin, time) {

    .validate_time_steps(time)

    q <- samc@p[-nrow(samc@p), -ncol(samc@p)]

    rdg <- as.vector(samc@p[-nrow(samc@p), ncol(samc@p)])

    time <- c(1, time)

    mort <- .sum_qpow_row(q, origin, time)

    mort <- lapply(mort, function(x){as.vector(x * rdg)})

    if (length(mort) == 1) {
      return(mort[[1]])
    } else {
      return(mort)
    }
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, dest, time) {

    .validate_time_steps(time)

    q <- samc@p[-nrow(samc@p), -ncol(samc@p)]

    rdg <- as.vector(samc@p[-nrow(samc@p), ncol(samc@p)])

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

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "numeric", time = "numeric"),
  function(samc, origin, dest, time) {

    mort <- mortality(samc, origin = origin, time = time)

    if (is.list(mort)){
      return(lapply(mort, "[", dest))
    } else if (is.vector(mort)) {
      return(mort[dest])
    } else {
      stop("Fatal error: This should not have been possible. Please submit a report with a fully reproducible and simplified example.")
    }
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "numeric"),
  function(samc, occ, time) {

    .validate_time_steps(time)

    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    Rdiag <- as.vector(samc@p[-nrow(samc@p), ncol(samc@p)])

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

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(mortality(samc, occ, time = time))
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc) {

    if (!samc@override)
      stop("This version of the mortality() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.")

    f <- visitation(samc)
    gc()
    rdg <- samc@p[-nrow(samc@p), ncol(samc@p)]
    r <- Matrix::sparseMatrix(i = 1:length(rdg),
                              j = 1:length(rdg),
                              x = rdg,
                              index1 = TRUE)

    # TODO f %*% r can be simplified to an elementwise multiplication of the matrix columns by the corresponding elements in the rdg vector. This might be helpful for memory allocations and performance.
    b <- f %*% r
    gc()
    return(b)
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "missing", time = "missing"),
  function(samc, origin) {

    rdg <- samc@p[-nrow(samc@p), ncol(samc@p)]

    vis <- visitation(samc, origin = origin)

    mort <- vis * rdg

    return(as.vector(mort))
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "missing"),
  function(samc, dest) {

    rdg <- samc@p[-nrow(samc@p), ncol(samc@p)]

    vis <- visitation(samc, dest = dest)

    mort <- vis * rdg[dest]

    return(as.vector(mort))
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "numeric", time = "missing"),
  function(samc, origin, dest) {

    rdg <- samc@p[-nrow(samc@p), ncol(samc@p)]

    vis <- visitation(samc, dest = dest)

    mort <- vis[origin] * rdg[dest]

    return(as.numeric(mort))
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {

    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    q = samc@p[-nrow(samc@p),-nrow(samc@p)]
    rdg <- as.vector(samc@p[-nrow(samc@p), ncol(samc@p)])

    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    pf <- .psif(q, pv)

    return(pf * rdg)
  })

#' @rdname mortality
setMethod(
  "mortality",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(mortality(samc, occ))
  })
