# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

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
#' This function is not intended to be used with other inputs such as `init` or `time`
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
        origin_is_matrix = is.matrix(origin)
        if (origin_is_matrix) {
            origin_df = as.data.frame(origin, stringsAsFactors = FALSE)
            if (ncol(origin_df) != 2) {
                stop("pairwise(): matrix origin must have exactly 2 columns.", call. = FALSE)
            }
            names(origin_df) = c("origin", "direction")
        } else if (is.atomic(origin) & is.vector(origin)) {
            origin_df = data.frame(origin = origin, stringsAsFactors = FALSE)
        } else {
            stop("pairwise(): only a vector or 2 column matrix supported for origin.", call. = FALSE)
        }

        dest_df = data.frame(dest = dest, stringsAsFactors = FALSE)

        df = merge(origin_df, dest_df, by = NULL) # Cartesian product: origin rows x dest values
        df = unique(df)
        rownames(df) = NULL

        origin_arg = if (origin_is_matrix) {
            as.matrix(df[, c("origin", "direction"), drop = FALSE])
        } else {
            df$origin
        }

        df$result = fun(samc, origin = origin_arg, dest = df$dest)
        df
    })

#' @rdname pairwise
setMethod(
    "pairwise",
    signature(fun = "function", samc = "samc", origin = "location", dest = "missing"),
    function(fun, samc, origin) {
        if (is.matrix(origin) | !is.atomic(origin)) {
            stop("pairwise(): when dest is excluded, origin must be a vector", call. = FALSE)
        }

        return(pairwise(fun, samc, origin, origin))
    })
