# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL


#' Rasterize input
#'
#' Convert input to a SpatRaster object
#'
#' This function is primarily used to convert an input matrix or raster to a SpatRaster
#' object. The main thing it is useful for is setting a standard extent and CRS for
#' converting matrices. It is used internally by the package to ensure consistent
#' results for the different data types for maps.
#'
#' When converting matrices, the extents are set to match the number of rows and
#' columns of the matrix. Pixels in the result are centered on whole number coordinates
#' with (1,1) corresponding to the bottom left pixel. The CRS is set to "local", which
#' treats it as Euclidean (Cartesian) plane with the units in meters.
#'
#' The main benefit will be for users that want an easy way to plot matrix data.
#' If the input type to the \code{\link{samc}} function is matrices, then the output
#' of \code{\link{map}} will also be matrices. Plotting these matrices can require
#' more work than simply using SpatRaster objects for \code{\link{samc}} and getting
#' SpatRaster results back from \code{\link{map}}.
#'
#' The raster and terra packages both also have a rasterize function that serves
#' a different purpose. If either of these packages are used directly, then the order of
#' package loading becomes very important because it will determine which version
#' of rasterize is used by default.
#'
#' @param x A matrix, RasterLayer, or SpatRaster
#
#' @return A SpatRaster object
#'
#' @example inst/examples/example.R
#'
#' @export
setGeneric(
  "rasterize",
  function(x) {
    standardGeneric("rasterize")
  })

#' @rdname rasterize
setMethod(
  "rasterize",
  signature(x = "matrix"),
  function(x) {
    terra::rast(x, extent = terra::ext(0.5, ncol(x) + 0.5, 0.5, nrow(x) + 0.5), crs = "local")
  })

#' @rdname rasterize
setMethod(
  "rasterize",
  signature(x = "RasterLayer"),
  function(x) {
    terra::rast(x)
  })

#' @rdname rasterize
setMethod(
  "rasterize",
  signature(x = "SpatRaster"),
  function(x) {
    x
  })
