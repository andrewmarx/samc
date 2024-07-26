# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.


#' Process location inputs
#'
#' Process location inputs
#'
#' @param samc A samc-class object
#' @param x A vector of integers or character names
#' @noRd
setGeneric(
  ".process_locations",
  function(samc, x) {
    standardGeneric(".process_locations")
  })

#' @noRd
setMethod(
  ".process_locations",
  signature(samc = "samc", x = "matrix"),
  function(samc, x) {
    if (samc@model$name == "CRW") {
      if (nrow(x) > 1) {
        stop("Multiple locations not supported yet. Matrix should only have 1 row/2 columns", call. = FALSE)
      }

      if (ncol(x) != 2) stop("Location should have 2 columns. The first for location and the second for direction.", call. = FALSE)

      .validate_locations(samc, x[1, 1])

      if (!(x[1, 2] %in% 1:8)) stop("Invalid direction. Must be a single integer from 1-8.", call. = FALSE)

      x = which(apply(samc@crw_map, 1, function(crw) return(all(crw == x))))

      if (length(x) != 1) stop("The combination of location and direction is not valid", call. = FALSE)

    } else {
      stop(paste("Invalid location input for model", samc@model$name), call. = FALSE)
    }
    return(x)
  })

#' @noRd
setMethod(
  ".process_locations",
  signature(samc = "samc", x = "numeric"),
  function(samc, x) {
    .validate_locations(samc, x)

    return(x)
  })

setMethod(
  ".process_locations",
  signature(samc = "samc", x = "character"),
  function(samc, x) {
    .validate_names(samc$names, x)

    return(.process_locations(samc, match(x, samc$names)))
  })


#' Process absorption inputs
#'
#' Process absorption inputs
#'
#' @param samc A samc-class object
#' @param x Absorption inputs
#' @noRd
setGeneric(
  ".process_abs_states",
  function(samc, x) {
    standardGeneric(".process_abs_states")
  })

#' @noRd
setMethod(
  ".process_abs_states",
  signature(samc = "samc", x = "SpatRaster"),
  function(samc, x) {

    x = rasterize(x)

    check(samc@map, x)

    if (terra::nlyr(x) == 0) {
      stop("Missing absorption data", call. = FALSE)
    }


    abs_mat = terra::values(x)
    abs_vec = as.vector(abs_mat[, 1])

    abs_minmax = terra::minmax(x)

    if (min(abs_minmax["min", ]) < 0 || max(abs_minmax["max", ]) > 1) {
      stop("Absorption values must be in range of 0-1", call. = FALSE)
    }

    excl = which(is.na(abs_vec))
    if (length(excl) > 0) {
      abs_mat = abs_mat[-excl, , drop = FALSE]
    }

    if (is.null(names(x))) colnames(abs_mat) = 1:ncol(abs_mat)

    if ("" %in% names(x)) stop("Mix of named and unnamed maps/layers", call. = FALSE)
    if (any(duplicated(names(x)))) stop("Duplicate names", call. = FALSE)

    colnames(abs_mat) = names(x)

    return(abs_mat)
  })

setMethod(
  ".process_abs_states",
  signature(samc = "samc", x = "RasterStack"),
  function(samc, x) {

    return(.process_abs_states(samc, terra::rast(x)))
  })

setMethod(
  ".process_abs_states",
  signature(samc = "samc", x = "RasterLayer"),
  function(samc, x) {

    return(.process_abs_states(samc, terra::rast(x)))
  })

setMethod(
  ".process_abs_states",
  signature(samc = "samc", x = "list"),
  function(samc, x) {

    # TODO: check against output field
    if (!all(sapply(x, is.matrix))) stop("List can only contain matrices. If using rasters, use raster::stack() or c() to combine terra SpatRasters instead.", call. = FALSE)

    x = lapply(x, rasterize)

    if(is.null(names(x))) names(x) = 1:length(x)

    return(.process_abs_states(samc, terra::rast(x)))
  })


#' Process initial state input
#'
#' Process initial state input
#'
#' @param samc A samc-class object
#' @param x initial state input
#' @noRd
setGeneric(
  ".process_init",
  function(samc, x) {
    standardGeneric(".process_init")
  })

# TODO: find a way to check the input type for `init` to the input type to samc()

#' @noRd
setMethod(
  ".process_init",
  signature(samc = "samc", x = "numeric"),
  function(samc, x) {
    if (any(!is.finite(x)) || any(x < 0)) stop("`init` input must only contain positive numeric values")

    if (length(x) != samc@nodes) stop("`init` input length does not match number of nodes")


    if (samc@model$name == "CRW") {
      x = sweep(samc@prob_mat[, terra::cells(samc@map)], 2, x, "*")
      dim(x) = NULL
      x = x[!is.na(x)]
    }

    return(x)
  })

#' @noRd
setMethod(
  ".process_init",
  signature(samc = "samc", x = "SpatRaster"),
  function(samc, x) {
    check(samc@map, x)

    pv = as.vector(terra::values(x))
    pv = pv[is.finite(pv)]

    return(.process_init(samc, pv))
  })

#' @noRd
setMethod(
  ".process_init",
  signature(samc = "samc", x = "RasterLayer"),
  function(samc, x) {
    return(.process_init(samc, rasterize(x)))
  })

setMethod(
  ".process_init",
  signature(samc = "samc", x = "matrix"),
  function(samc, x) {
    return(.process_init(samc, rasterize(x)))
  })


#' Map Location
#'
#' Map a location
#'
#' @param x A list
#' @noRd
.map_location = function(samc, x) {
  if (samc@source == "transition") {
    vec = numeric(samc@nodes)
    vec[x] = 1
    names(vec) = samc@names

    return(vec)
  } else {
    df = data.frame(cell = terra::cells(samc@map),
                    vec = numeric(samc@nodes))
    df$vec[x] = 1
  }

  return(.build_map(samc, df))
}


#' Summarize CRW
#'
#' Summarize CRW
#'
#' @param samc samc model
#' @noRd
.summarize_crw = function(samc, vec, fun) {
  if (length(vec) != nrow(samc@crw_map))
    stop("The length of the vector does not match the number of non-NA cells in the landscape data", call. = FALSE)

  aggregate(vec ~ samc@crw_map[, 1], FUN = fun)$vec
}
