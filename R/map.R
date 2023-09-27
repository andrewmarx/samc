# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include check.R samc-class.R
NULL


#' Map vector data
#'
#' Map vector data to a RasterLayer
#'
#' This is a convenience function to ensure that vector data is properly mapped
#' back to the original landscape data. The reason this is needed is that the package
#' supports matrices, RasterLayers, and SpatRasters, which can differ in the order
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
#' @return A matrix, RasterLayer, or SpatRaster object. The returned type will match
#' the type used to create the samc object.
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
  signature(samc = "samc", vec = "numeric"),
  function(samc, vec) {
    # TODO make work for transition matrices
    if (samc@source == "transition") stop("This function cannot be used for samc objects created from transition matrices", call. = FALSE)

    if (length(vec) != length(terra::cells(samc@map)))
      stop("The length of the vector does not match the number of non-NA cells in the landscape data", call. = FALSE)

    df = data.frame(cell = terra::cells(samc@map),
                    vec = vec)

    .build_map(samc, df)
  })

#' @rdname map
setMethod(
  "map",
  signature(samc = "samc", vec = "list"),
  function(samc, vec){

    lapply(vec, function(x){
      map(samc, x)
    })
  })


#' Build map
#'
#' Internal function to map a df to various objects
#'
#' @param samc samc
#' @param df df
#' @noRd
.build_map = function(samc, df) {
  ras_base = as.numeric(samc@map)

  ras_list = lapply(sort(colnames(df)[-1]), function (x) {
    ras = ras_base
    ras[terra::cells(ras)] = df[x]

    return(ras)
  })

  ras = terra::rast(ras_list)

  if (samc@source == "SpatRaster") {
    return(ras)
  } else if (samc@source == "RasterLayer") {
    if (terra::nlyr(ras) > 1) {
      return(raster::stack(ras))
    } else {
      return(raster::raster(ras))
    }
  } else if (samc@source == "matrix") {
    if (terra::nlyr(ras) > 1) {
      return(lapply(ras, function(x) {
        as.matrix(x, wide = TRUE)
      }))
    } else {
      return(as.matrix(ras, wide = TRUE))
    }
  } else {
    stop("An unexpected error occurred. Please report as a bug with a reproducible example", call. = FALSE)
  }
}
