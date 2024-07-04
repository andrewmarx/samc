# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL


#' Check landscape data
#'
#' Check that landscape inputs have valid values and matching properties.
#'
#' This function is used to ensure that inputs (resistance, absorption, fidelity,
#' and occupancy) have valid values and the same properties. This includes
#' checking the CRS (if using raster inputs), dimensions, and locations of
#' cells with NA data. It can be used to directly compare two matrices or two
#' rasters, or it can be used to check a \code{\link{samc-class}} object
#' against a matrix or raster.
#'
#' It can also be used to check a numeric vector against a \code{\link{samc-class}} object
#' created from a P matrix. In this case, the length of the vector must be equal to
#' the number of transient states. If the transient states are named, the vector
#' must contain the same names.
#'
#' The function returns \code{TRUE} if the inputs have matching properties. Otherwise,
#' it will stop execution and print an error message with details about the
#' difference between the two inputs.
#'
#' Note that the package assumes the different landscape inputs will be the same
#' type, either matrices or RasterLayers. Mixing RasterLayer data and matrix
#' data is not supported.
#'
#' @param a A \code{\link{samc-class}}, \code{\link{matrix}}, or \code{\link[raster]{RasterLayer-class}} object
#' @param b A \code{\link{matrix}} or \code{\link[raster]{RasterLayer-class}} object
#
#' @return See \emph{Details} section.
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "check",
  function(a, b) {
    standardGeneric("check")
  })


# TODO merge
#' @rdname check
setMethod(
  "check",
  signature(a = "Raster", b = "missing"),
  function(a){

    n <- raster::nlayers(a)

    if (n == 0) {
      stop("No raster layers found", call. = FALSE)
    }

    r1 <- a[[1]]

    if (sum(is.infinite(r1[]), na.rm = TRUE) > 0) {
      stop("Data contains Inf or -Inf element", call. = FALSE)
    } else if (sum(is.nan(r1[]), na.rm = TRUE) > 0) {
      stop("Data contains NaN elements", call. = FALSE)
    }

    if (n > 1) {
      r1[] <- is.finite(r1[])

      for (i in 2:n) {
        r2 <- a[[i]]

        if (sum(is.infinite(r2[]), na.rm = TRUE) > 0) {
          stop("Data contains Inf or -Inf element", call. = FALSE)
        } else if (sum(is.nan(r2[]), na.rm = TRUE) > 0) {
          stop("Data contains NaN elements", call. = FALSE)
        }

        r2[] <- is.finite(r2[])

        tryCatch(
          {
            raster::compareRaster(r1, r2, values = TRUE)
          },
          error = function(e) {
            if(grepl("not all objects have the same values", e$message)) {
              msg = "NA mismatch"
            } else {
              msg = e$message
            }
            stop(msg, " in input data", call. = FALSE)
          }
        )
      }
    }

    return(TRUE)
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "SpatRaster", b = "missing"),
  function(a){

    n <- terra::nlyr(a)

    if (n == 0) {
      stop("No raster layers found", call. = FALSE)
    }

    inf_counts = numeric(n)
    nan_counts = numeric(n)

    for (r in 1:terra::nrow(a)) {
      data = terra::values(a, row = r, nrows = 1)

      if (any(is.infinite(data))) stop("Data contains Inf or -Inf", call. = FALSE)
      if (any(is.nan(data))) stop("Data contains Inf or -Inf", call. = FALSE)

      data = rowSums(is.finite(data))
      if (any(data > 0 & data < n)) stop("NA mismatch in input data", call. = FALSE)
    }


    return(TRUE)
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "matrix", b = "missing"),
  function(a){
    a <- rasterize(a)

    check(a)
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "SpatRaster", b = "SpatRaster"),
  function(a, b){
    check(c(a, b)) # TODO make CRS warning an error?
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "Raster", b = "Raster"),
  function(a, b){
    check(raster::stack(a, b))
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "matrix", b = "matrix"),
  function(a, b){
    a <- rasterize(a)
    b <- rasterize(b)

    check(a, b)
  })


# TODO: reimplement source checks to ensure inputs match original data type used to create samc object

#' @rdname check
setMethod(
  "check",
  signature(a = "samc", b = "Raster"),
  function(a, b){
    check(a@map, b)
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "samc", b = "SpatRaster"),
  function(a, b){
    check(a@map, b)
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "samc", b = "matrix"),
  function(a, b){
    check(a@map, rasterize(b))
  })


#' @rdname check
setMethod(
  "check",
  signature(a = "samc", b = "numeric"),
  function(a, b){
    if (a@source != "transition") stop("Numeric vector input only valid for samc objects created from P matrix", call. = FALSE)

    if (!isTRUE(all.equal(names(b), a@names))) stop("Names of the vector must match the names of the transient states in the P matrix", call. = FALSE)

    if (any(!is.finite(b)) || any(b < 0)) stop("Input must only contain positive numeric values", call. = FALSE)

    if (length(b) != length(a@data@t_abs)) stop("Input length does not match number of transient states", call. = FALSE)
  })
