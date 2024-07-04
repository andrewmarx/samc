# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R visitation.R
NULL


#' Calculate survival metrics
#'
#' Calculates the expected time to absorption
#'
#' \eqn{z=(I-Q)^{-1}{\cdot}1=F{\cdot}1}
#' \itemize{
#'   \item \strong{survival(samc)}
#'
#' The result is a vector \eqn{\mathbf{v}} where \eqn{\mathbf{v}_i} is the expected
#' time to absorption if starting at transient state \eqn{\mathit{i}}.
#'
#' If the samc-class object was created using matrix or RasterLayer maps, then
#' vector \eqn{\mathbf{v}} can be mapped to a RasterLayer using the
#' \code{\link{map}} function.
#' }
#'
#' \eqn{\psi^Tz}
#' \itemize{
#'   \item \strong{survival(samc, init)}
#'
#' The result is a numeric that is the expected time to absorption given an initial
#' state \eqn{\psi}.
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-init
#'
#' @return See Details
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "survival",
  function(samc, init, origin) {
    standardGeneric("survival")
  })

# survival(samc) ----
#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", init = "missing", origin = "missing"),
  function(samc) {
    .disable_conv(samc)

    # TODO If multiple dest becomes supported, this can be replaced with visitation(dest)
    if (samc@solver == "iter") {
      z = .f1_iter(samc@data@f)
    } else {
      z = .f1(samc@data@f, samc@.cache$sc)
    }

    if (samc@model$name == "CRW") {
      vec = as.vector(samc@prob_mat)
      vec = vec[!is.na(vec)]

      z = vec * z
      z = .summarize_crw(samc, z, sum)
    }

    return(z)
  })

# survival(samc, origin) ----
#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", init = "missing", origin = "location"),
  function(samc, origin) {
    if (is(origin, "matrix")) {
      if (nrow(origin) > 1) stop("Only a single origin is supported for CRW", call. = FALSE)
    } else {
      if (length(origin) != 1)
        stop("origin can only contain a single value for this version of the function", call. = FALSE)
    }

    origin = .process_locations(samc, origin)
    init = .map_location(samc, origin)

    return(survival(samc, init = init))
  })

# survival(samc, init) ----
#' @rdname survival
setMethod(
  "survival",
  signature(samc = "samc", init = "ANY", origin = "missing"),
  function(samc, init) {
    sum(visitation(samc, init))
  })
