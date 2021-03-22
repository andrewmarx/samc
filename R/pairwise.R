# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details..

#' @include samc-class.R location-class.R
NULL

#' Pairwise analyses
#'
#' Analysis for pairwise combinations locations
#'
#' When providing vector inputs for the `origin` and `dest` parameters to analytical
#' functions, the package assumes that users are providing pairs of `origin` and
#' `dest`. That is, `origin[1]` is paired with `dest[1]`, `origin[2]` is paired
#' `dest[2]`, etc. Another way to think about it is that these two vector inputs
#' can be treated as columns in the same dataframe. The result of the analytical
#' function then is a vector of the same length as the input. This behavior works
#' for any situation, so it is the default for the package.
#'
#' However, some users may wish to run an analytical function for all the pairwise
#' combinations of the values in the input vectors. That is, `origin[1]` is paired
#' with `dest[1]`,`dest[2]`, `dest[3]`, etc, before moving on to the next elements
#' in `origin`. This approach has the advantage of potentially reducing the amount
#' of code needed for an analysis, and the results can be represented as a pairwise
#' matrix, but it is not suitable for all situations. To enable this second approach
#' more easily, the `pairwise()` function runs all the combinations of the `origin`
#' and `dest` parameters for an analytical function and returns the results in a
#' 'long' format data.frame. This data.frame can then be reshaped into a pairwise
#' matrix or 'wide' format data.frame using tools like the reshape2 or tidyr packages.
#'
#' This function is not intended to be used with other inputs such as `occ` or `time`
#'
#' @param fun A samc analytical function with signature fun(samc, origin, dest)
#' @param samc A \code{\link{samc-class}} object
#' @param origin A vector of locations
#' @param dest A vector of locations. Can be excluded to reuse the origin parameter
#'
#' @return A 'long' format data.frame
#'
#' @example inst/examples/pairwise.R
#'
#' @export

setGeneric(
  "pairwise",
  function(fun, samc, origin, dest) {
    standardGeneric("pairwise")
  })

#' @rdname pairwise
setMethod(
  "pairwise",
  signature(fun = "function", samc = "samc", origin = "location", dest = "location"),
  function(fun, samc, origin, dest) {
    # Create all possible pairs
    df <- expand.grid(origin = origin, dest = dest,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

    # Remove duplicates
    df <- unique(df)

    # Reset rownames
    rownames(df) <- NULL

    # Get results
    df$result <- fun(samc, origin = df$origin, dest = df$dest)

    return(df)
  })

#' @rdname pairwise
setMethod(
  "pairwise",
  signature(fun = "function", samc = "samc", origin = "location", dest = "missing"),
  function(fun, samc, origin) {
    return(pairwise(fun, samc, origin, origin))
  })
