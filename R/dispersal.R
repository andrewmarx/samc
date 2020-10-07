# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate dispersal metrics
#'
#' Calculates the probability of individuals visiting locations
#'
#'
#' \eqn{\tilde{D}_{jt}=(\sum_{n=0}^{t-1}\tilde{Q}^n)\tilde{q}_j}
#' \itemize{
#'   \item \strong{dispersal(samc, dest, time)}
#'
#' The result is a vector (single time step) or a list of vectors (multiple
#' time steps) where each element corresponds to a cell in the landscape,
#' and can be mapped back to the landscape using the \code{\link{map}} function.
#' Element \emph{k} is the probability of ever visiting a given destination,
#' if starting at any other location, within \emph{t} or fewer time steps.
#' }
#'
#' \eqn{\psi^T\tilde{D}_{jt}}
#' \itemize{
#'   \item \strong{dispersal(samc, occ, dest, time)}
#'
#' The result is a numeric (single time step) or a list of numerics (multiple
#' time steps) that is the unconditional probability of visiting a
#' given destination within \emph{t} or fewer time steps.
#' }
#'
#' \eqn{D=(F-I)diag(F)^{-1}}
#' \itemize{
#'   \item \strong{dispersal(samc)}
#'
#' The result is a matrix where element (\emph{i},\emph{j}) is the probability
#' that location \emph{j} is visited when starting in location \emph{i}.
#'
#' The returned matrix will always be dense and cannot be optimized. Must enable
#' override to use.
#'
#'   \item \strong{dispersal(samc, origin)}
#'
#' This function has not been optimized yet, and will not run.
#'
#'   \item \strong{dispersal(samc, dest)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element \emph{i} is the probability that the
#' destination is visited when starting in location \emph{i}.
#'
#'   \item \strong{dispersal(samc, origin, dest)}
#'
#' The result is a numeric value that is the probability that an individual
#' starting at the origin visits the destination.
#' }
#'
#' \eqn{\psi^TD}
#' \itemize{
#'   \item \strong{dispersal(samc, occ)}
#'
#' The result is a vector where each element corresponds to a cell in the
#' landscape, and can be mapped back to the landscape using the
#' \code{\link{map}} function. Element \emph{j} is the unconditional
#' probability distribution of ever visiting location \emph{j}, regardless of
#' the initial location.
#'
#'   \item \strong{dispersal(samc, occ, dest)}
#'
#' The result is a numeric value that is the unconditional probability
#' distribution of ever visiting a given destination, regardless of the initial
#' location.
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
#' @return A matrix, a vector, a list of vectors, or a numeric
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "dispersal",
  function(samc, occ, origin, dest, time) {
    standardGeneric("dispersal")
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, dest, time) {

    .validate_time_steps(time)

    if (dest %% 1 != 0 || dest < 1 || dest > sum(samc@map[], na.rm = TRUE))
      stop("dest must be an integer that refers to a cell in the landscape")

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    qv <- q[, dest]
    qv <- qv[-dest]
    q <- q[-dest, -dest]

    q2 = q
    q2@x <- -q2@x
    Matrix::diag(q2) <- Matrix::diag(q2) + 1

    time <- c(0, time)
    res <- .sum_qn_q(q, q2, qv, time)

    res <- lapply(res, as.vector)

    if (length(res) == 1) {
      return(res[[1]])
    } else {
      return(res)
    }
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, occ, dest, time) {

    check(samc, occ)

    d <- dispersal(samc, dest = dest, time = time)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]
    pv <- pv[-dest]

    if (is.list(d)) {
      return(lapply(d, FUN = function(x){as.numeric(pv %*% x)}))
    } else {
      return(as.numeric(pv %*% d))
    }
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "numeric", time = "numeric"),
  function(samc, occ, dest, time) {

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(dispersal(samc, occ, dest = dest, time = time))
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "missing", time = "missing"),
  function(samc) {

    if (!samc@override)
      stop("This version of the dispersal() method produces a large dense matrix.\nIn order to run it, create the samc object with the override parameter set to TRUE.")

    f <- visitation(samc)
    gc()
    fdg <- 1 / Matrix::diag(f)
    fdg_mat <- Matrix::sparseMatrix(i = 1:length(fdg),
                                    j = 1:length(fdg),
                                    x = fdg,
                                    index1 = TRUE)
    # TODO Check if 'diag(f) <- ' resulting in extra memory allocation.
    Matrix::diag(f) <- Matrix::diag(f) - 1
    gc()
    d_mat <- f %*% fdg_mat

    return(d_mat)
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "missing", time = "missing"),
  function(samc, origin) {
    stop("A suitably optimized version of this function has not been identified (yet). As a workaround, consider calculation destination columns instead")
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "missing"),
  function(samc, dest) {
    f_col <- visitation(samc, dest = dest)
    fjj <- f_col[dest]
    f_col[dest] <- f_col[dest] - 1

    d_vec <- f_col/fjj

    return(d_vec)
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "missing", origin = "numeric", dest = "numeric", time = "missing"),
  function(samc, origin, dest) {
    d <- dispersal(samc, dest = dest)

    return(d[origin])
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {

    check(samc, occ)

    q <- samc@p[-nrow(samc@p), -nrow(samc@p)]
    q@x <- -q@x
    Matrix::diag(q) <- Matrix::diag(q) + 1

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    disp <- .psid_long(q, pv)

    return(as.vector(disp))
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "missing", time = "missing"),
  function(samc, occ) {

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(dispersal(samc, occ))
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "RasterLayer", origin = "missing", dest = "numeric", time = "missing"),
  function(samc, occ, dest) {

    check(samc, occ)

    pv <- as.vector(occ)
    pv <- pv[is.finite(pv)]

    dj <- dispersal(samc, dest = dest)

    return(as.numeric(pv %*% dj))
  })

#' @rdname dispersal
setMethod(
  "dispersal",
  signature(samc = "samc", occ = "matrix", origin = "missing", dest = "numeric", time = "missing"),
  function(samc, occ, dest) {

    occ <- raster::raster(occ, xmn = 0.5, xmx = ncol(occ) + 0.5, ymn = 0.5, ymx = nrow(occ) + 0.5)

    return(dispersal(samc, occ, dest = dest))
  })
