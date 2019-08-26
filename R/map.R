# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include check.R samc-class.R
NULL


#' Map vector data
#'
#' Map vector data to a RasterLayer
#'
#' This is a convenience function to ensure that vector data is properly mapped
#' back to the original landscape data. The reason this is needed is that the
#' package supports both matrices and RasterLayers, which differ in the order
#' that data is read and written (R matrices are column-major order, whereas the
#' raster package uses row-major order). Internally, the package uses only a
#' single order, regardless of the original data. This can cause issues with
#' mapping vector results if care is not taken, and this function is provided to
#' simplify the process. It also correctly maps results for landscape data that
#' has NA cells, which are another potential source of error if not careful.
#'
#' The only requirement of the \code{vec} input is that the number of elements
#' in it matches the number of non-NA cells in the landscape data that was used
#' to create the samc object.
#'
#' @param samc Spatial absorbing Markov chain object. This should be output from the samc() function.
#' @param vec Vector data to fill into the map.
#' @return A RasterLayer object
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "map",
  function(samc, vec) {
    standardGeneric("map")
  })

#' @rdname map
setMethod(
  "map",
  signature(samc = "samc", vec = "vector"),
  function(samc, vec){

    if (length(vec) != sum(samc@map[], na.rm = TRUE))
      stop("The length of the vector does not match the number of non-NA cells in the landscape data")

    ras <- samc@map

    ras[ras[]] <- vec
    ras[!samc@map[]] <- NA

    return(ras)
  })
