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
#'   \item \strong{cond_passage(samc, origin, dest)}
#'
#' The result is a numeric value representing the mean number of steps before
#' absorption starting from a given origin conditional on absorption into \emph{j}.
#' }
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

    t <- .cond_t(Qj, qj)

    return(as.numeric(t))
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
      adj_origin <- origin
      adj_origin[origin > d] <- adj_origin[origin > d] - 1
      result[dest == d] <- t[adj_origin[dest == d]]
    }

    result[dest == origin] <- NA

    return(result)
  })
