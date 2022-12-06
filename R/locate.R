# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details..

#' @include samc-class.R
NULL


#' Get cell numbers
#'
#' Get cell numbers from raster data
#'
#' This function is used to get cell numbers from raster data. The numbers used
#' for origin and destination values in many samc metrics refer to column/row
#' numbers of the P matrix. For a P matrix derived from raster data, these numbers
#' would normally line up with the cell numbers of the raster, but this is not
#' always true. This is the case when the raster contains NA data; the cells
#' associated with this data are excluded from the P matrix. This causes issues
#' trying to determine the cell numbers that should be used in analyses.
#'
#' The \code{\link{locate}} function operates more-or-less like the
#' \code{\link[raster]{cellFromXY}} function in the raster package, but unlike
#' \code{\link[raster]{cellFromXY}}, locate properly accounts for NA cells
#' in identifying cell numbers from coordinate data.
#'
#' This function can also be used if the samc object was created from matrix inputs
#' for the resistance, absorption, and fidelity parameters. In this case, the
#' values in the xy coordinate parameter can be column-row values with the caveat
#' that (1,1) is the bottom left corner.
#'
#' The xy parameter can also be excluded. In this case, the function returns a
#' raster where the values of the cells contains the cell number.
#'
#' Internally, this function relies on the \code{\link[raster]{extract}} function
#' from the raster package, and any valid input for the y argument of that function
#' is valid here.
#'
#' @param samc A \code{\link{samc-class}} object
#' @param xy Any valid input to the y argument of the \code{\link[raster]{extract}} function in the raster package.
#
#' @return A RasterLayer or a vector
#'
#' @example inst/examples/locate.R
#'
#' @export

setGeneric(
  "locate",
  function(samc, xy) {
    standardGeneric("locate")
  })

#' @rdname locate
setMethod(
  "locate",
  signature(samc = "samc", xy = "missing"),
  function(samc){
    if (samc@source == "transition") {
      stop("This function cannot be used when the samc-object was created from a transition matrix", call. = FALSE)
    } else if (samc@source == "SpatRaster") {
      return(samc@map)
    } else if (samc@source == "RasterLayer") {
      return(raster::raster(samc@map))
    } else if (samc@source == "matrix") {
      return(as.matrix(samc@map, wide = TRUE))
    } else {
      stop("An unexpected issue occurred. Please report a bug with a reproducible example", call. = FALSE)
    }
  })

#' @rdname locate
setMethod(
  "locate",
  signature(samc = "samc", xy = "ANY"),
  function(samc, xy){
    if (samc@source == "transition") stop("This function cannot be used when the samc-object was created from a transition matrix", call. = FALSE)

    result <- terra::extract(samc@map, xy)[, 2]

    if (anyNA(result)) stop("One or more coordinates do not correspond to non-NA cells.", call. = FALSE)

    return(result)
  })
