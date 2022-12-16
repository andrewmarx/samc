# Copyright (c) 2022 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.


#' map class
#'
#' Union class for map inputs
#'
#' The map class is a union class of the "matrix", "RasterLayer", and "SpatRaster"
#' classes. Users generally do not need to worry about it except to know that any
#' method parameter with "map" as the type can have one of these three classes
#' provided as input.
#'
#' @export map

setClassUnion("map", c("matrix", "RasterLayer", "SpatRaster"))
