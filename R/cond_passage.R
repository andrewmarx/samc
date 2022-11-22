# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R location-class.R
NULL

#' Conditional Mean First Passage Time
#'
#' Calculate the mean number of steps to first passage
#'
#'
#' \eqn{\tilde{t}=\tilde{B}_j^{-1}\tilde{F}\tilde{B}_j{\cdot}1}
#' \itemize{
#'   \item \strong{cond_passage(samc, dest)}
#'
#' The result is a vector where each element corresponds to a cell in the landscape,
#' and can be mapped back to the landscape using the \code{\link{map}} function.
#' Element \emph{i} is the mean number of steps before absorption starting from
#' location \emph{i} conditional on absorption into \emph{j}
#'
#' Note that mathematically, the formula actually does not return a value for when
#' \emph{i} is equal to \emph{j}. This leads to a situation where the resultant vector
#' is actually one element short and the index for some of the elements may be shifted.
#' The \strong{cond_passage()} function fills inserts a \code{0} value for vector indices
#' corresponding to \emph{i == j}. This corrects the final result so that vector indices
#' work as expected, and allows the vector to be properly used in the \code{\link{map}} function.
#'
#'   \item \strong{cond_passage(samc, origin, dest)}
#'
#' The result is a numeric value representing the mean number of steps before
#' absorption starting from a given origin conditional on absorption into \emph{j}.
#'
#' As described above, mathematically the formula does not return a result for when
#' the \code{origin} and \code{dest} inputs are equal, so the function simply returns a \code{0}
#' in this case.
#' }
#'
#' \strong{WARNING}: This function is not compatible when used with data where there
#' are states with total absorption present. When present, states representing
#' total absorption leads to unsolvable linear equations. The only exception to this
#' is when there is a single total absorption state that corresponds to input to
#' the \code{dest} parameter. In this case, the total absorption is effectively
#' ignored when the linear equations are solved.
#'
#' \strong{WARNING}: This function will crash when used with data representing
#' a disconnected graph. This includes, for example, isolated pixels or islands
#' in raster data. This is a result of the transition matrix for disconnected
#' graphs leading to some equations being unsolvable. Different options
#' are being explored for how to best identify these situations in data and
#' handle them accordingly.
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-origin
#' @template param-dest
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "cond_passage",
  function(samc, origin, dest) {
    standardGeneric("cond_passage")
  })

# cond_passage(samc, dest) ----
#' @rdname cond_passage
setMethod(
  "cond_passage",
  signature(samc = "samc", origin = "missing", dest = "location"),
  function(samc, dest) {
    if (samc@clumps == -1)
      warning("Unknown number of clumps in data. If the function crashes, it may be due to the transition matrix being discontinuous.", call. = FALSE)

    if (samc@clumps > 1)
      stop("This function cannot be used with discontinuous data", call. = FALSE)

    if (length(dest) != 1)
      stop("dest must be a single location that refers to a cell in the landscape", call. = FALSE)

    dest <- .process_locations(samc, dest)

    Q <- samc$q_matrix
    qj <- Q[-dest, dest]
    Qj <- Q[-dest, -dest]

    Qj@x <- -Qj@x
    Matrix::diag(Qj) <- Matrix::diag(Qj) + 1

    if (samc@solver == "iter") {
      t <- as.numeric(.cond_t_iter(Qj, qj))
    } else {
      t <- as.numeric(.cond_t(Qj, qj))
    }

    # insert 0 element back into vector so output length is same original data
    final <- 1:(length(t) + 1)
    names(final) <- samc$names

    j <- 1
    for (i in 1:length(final)) {
      if (final[dest] == i) {
        final[dest] <- 0
        next
      }
      final[i] <- t[j]
      j <- j + 1
    }

    return(final)
  })

# cond_passage(samc, origin, dest) ----
#' @rdname cond_passage
setMethod(
  "cond_passage",
  signature(samc = "samc", origin = "location", dest = "location"),
  function(samc, origin, dest) {
    if(length(origin) != length(dest))
      stop("The 'origin' and 'dest' parameters must have the same number of values", call. = FALSE)

    origin <- .process_locations(samc, origin)
    dest <- .process_locations(samc, dest)

    result <- vector(mode = "numeric", length = length(origin))

    unique_dest <- unique(dest)

    for (d in unique_dest) {
      t <- cond_passage(samc, dest = d)
#      adj_origin <- origin
#      adj_origin[origin > d] <- adj_origin[origin > d] - 1
      result[dest == d] <- t[origin[dest == d]]
    }

    return(result)
  })
