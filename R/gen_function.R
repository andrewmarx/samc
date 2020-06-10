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
#'   \item \strong{gen_function(samc)}
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
  "gen_function",
  function(samc) {
    standardGeneric("gen_function")
  })

#' @rdname gen_function
setMethod(
  "gen_function",
  signature(samc = "samc", occ = "missing", origin = "missing", dest = "numeric", time = "missing"),
  function(samc, dest) {

  })
