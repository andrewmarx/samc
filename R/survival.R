# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

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

    if (samc@solver == "iter") {
      z = .f1_iter(samc@data@f)
    } else {
      z = .f1(samc@data@f, samc@.cache$sc)
    }

    z = as.vector(z)

    if (samc@model$name == "CRW") {
      z = .summarize_crw(samc, z, sum) # TODO figure out why this produces slightly incorrect results
      warning("survival() CRW results can be incorrect by up to 0.5%") # TODO fix issue and remove warning
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
    if (samc@solver %in% c("direct", "iter")) {
      check(samc, init)

      pv <- .process_init(samc, init)

      sv <- survival(samc)

      if (samc@model$name == "CRW") {
        pv = .summarize_crw(samc, pv, sum)
        warning("survival() CRW results can be incorrect by up to 0.5%") # TODO fix issue and remove warning
      }

      surv <- pv %*% sv

      return(as.numeric(surv))
    } else if (samc@solver == "conv") {

      res = visitation(samc, init)

      return(sum(res))
    } else {
      stop("Invalid method attribute in samc object.")
    }
  })
