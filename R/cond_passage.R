# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL

#' Short title
#'
#' Short description
#'
#'
#' \eqn{\bar{t}_p=diag(\tilde{B})^{-1}\tilde{F}diag(\tilde{B}){\cdot}1}
#' \itemize{
#'   \item \strong{cond_passage(samc, origin, dest)}
#'
#' Long description
#' }
#'
#' @template section-perf
#'
#' @template param-samc
#' @template param-dest
#'
#' @return Result description
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "cond_passage",
  function(samc, origin, dest) {
    standardGeneric("cond_passage")
  })

#' @rdname cond_passage
setMethod(
  "cond_passage",
  signature(samc = "samc", origin = "missing", dest = "numeric"),
  function(samc, dest) {

  })
