# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.
#
# This file is for internal functions. They are subject to change and should not
# be used by users.
#

#' Transition function
#'
#' An internal function for creating transition matrices
#'
#' @param data A SpatRaster
#' @noRd
.transition <- function(resistance, absorption, fidelity, fun, dir, sym = TRUE) {
  if (class(fun) == "character" || !(dir %in% c(4, 8))) {
    stop("gdistance's named funtion options not supported")
    #return(gdistance::transition(data, fun, dir, sym))
  }

  #data_crs = terra::crs(resistance)
  #data_cells = terra::ncell(resistance)
  #data_rows = terra::nrow(resistance)
  #data_cols = terra::ncol(resistance)

  cell_nums = terra::cells(resistance)
  adj = terra::adjacent(resistance, cells=cell_nums, pairs=TRUE, directions=dir)
  adj = adj[adj[, 2] %in% cell_nums, ]

  adj_length = nrow(adj)

  if(sym) adj = adj[adj[, 1] < adj[, 2], ]

  coords <- cbind(terra::xyFromCell(resistance, adj[, 1]),
                  terra::xyFromCell(resistance, adj[, 2]))


  dist <- terra::distance(coords[,1:2],
                          coords[,3:4],
                          lonlat = terra::is.lonlat(resistance),
                          pairwise = TRUE)


  i = adj[,1]
  j = adj[,2]

  sums = numeric(max(adj))

  res = as.vector(terra::values(resistance))

  adj[] = res[adj]
  rm(res);gc()

  transition.values = numeric(adj_length)

  if (sym) {
    adj_length = nrow(adj)
    for (k in 1:adj_length) {
      transition.values[c(k, k + adj_length)] = fun(adj[k, ]) / dist[k]

      sums[i[k]] = sums[i[k]] + transition.values[k]
      sums[j[k]] = sums[j[k]] + transition.values[k]
    }

    tmp = i
    i = c(i, j)
    j = c(j, tmp)
  } else {
    for (k in 1:nrow(adj)) {
      transition.values[k] = fun(adj[k, ]) / dist[k]

      sums[i[k]] = sums[i[k]] + transition.values[k]
    }
  }

  rm(adj);rm(dist);gc()

  tmp = as.vector(terra::values(1 - absorption - fidelity))

  transition.values = tmp[i] * transition.values / sums[i]

  rm(tmp, sums); gc()

  if(!all(transition.values >= 0)){
    warning("transition function gives negative values")
  }

  i = c(i, cell_nums)
  j = c(j, cell_nums)
  transition.values = c(transition.values, as.vector(terra::values(fidelity))[cell_nums])


  # Adjust for NAs
  i = match(i, cell_nums)
  j = match(j, cell_nums)


  transitionMatrix = Matrix::sparseMatrix(i = i, j = j, x = transition.values)

  return(transitionMatrix)
}


#' Validate time steps
#'
#' Performs several checks to make sure a vector of time steps is valid
#'
#' @param x A vector object to be validated as time steps
#' @noRd
.validate_time_steps <- function(x) {
  if (!is.numeric(x))
    stop("The time argument must be a positive integer or a vector of positive integers", call. = FALSE)

  if (sum(is.na(x)) > 0)
    stop("NA values are not allowed in the time argument", call. = FALSE)

  if (any(x %% 1 != 0))
    stop("Decimal values are not allowed in the time argument", call. = FALSE)

  if (any(x < 1))
    stop("All time steps must be positive (greater than 0)", call. = FALSE)

  if (is.unsorted(x))
    stop("The provided time steps must be in ascending order.", call. = FALSE)

  if (sum(duplicated(x) > 0))
    stop("Duplicate time steps are not allowed in the time argument", call. = FALSE)

  if (any(x > 10000))
    stop("Due to how the short-term metrics are calculated and the way that
  decimal numbers are handled by computers, numerical issues related to
  precision arise when a time step value is too high. Currently, a hard
  limit of 10000 time steps is enforced to encourage users to more
  seriously consider how many time steps are relevant to their use case.
  For example, if a single time step represents 1 day, then the current
  limit represents 24.7 years. There is flexibility to increase the limit
  if a justification can be made for it, but it's far more likely that
  users will generally want far fewer time steps for ecologically relevant
  results and to avoid the cummulative precision issues.", call. = FALSE)
}


#' Validate location vectors
#'
#' Performs several checks to make sure a vector locations is valid
#'
#' @param samc samc-class object
#' @param x A vector object to be validated as locations
#' @noRd
.validate_locations <- function(samc, x) {
  if (!is.numeric(x))
    stop("Locations must be a positive integer or a vector of positive integers", call. = FALSE)

  if (sum(is.na(x)) > 0)
    stop("NA values are not valid locations", call. = FALSE)

  if (any(x %% 1 != 0))
    stop("Decimal values are not valid locations", call. = FALSE)

  if (any(x < 1))
    stop("All location values must be positive (greater than 0)", call. = FALSE)

  if (any(x > nrow(samc$q_matrix)))
    stop("Location values cannot exceed the number of nodes in the landscape", call. = FALSE)
}


#' Validate location names
#'
#' Performs several checks to make sure a vector of names is valid
#'
#' @param vec A vector of location names
#' @param x A vector object to be validated as names
#' @noRd
.validate_names <- function(vec, x) {
  invalid_names <- x[!(x %in% vec)]

  if (length(invalid_names > 0)){
    print(vec)
    print(x)
    stop(paste("\nInvalid location name:", invalid_names), call. = FALSE)
  }
}


#' Rasterize matrices
#'
#' Convert a matrix to a SpatRaster. Ensures consistency of conversion throughout the package
#'
#' @param x A matrix
#' @noRd
setGeneric(
  ".rasterize",
  function(x) {
    standardGeneric(".rasterize")
  })

#' @noRd
setMethod(
  ".rasterize",
  signature(x = "matrix"),
  function(x) {
    terra::rast(x, extent = terra::ext(0.5, ncol(x) + 0.5, 0.5, nrow(x) + 0.5), crs = "local")
  })

#' @noRd
setMethod(
  ".rasterize",
  signature(x = "RasterLayer"),
  function(x) {
    terra::rast(x)
  })

#' @noRd
setMethod(
  ".rasterize",
  signature(x = "SpatRaster"),
  function(x) {
    x
  })



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

    return(match(x, samc$names))
  })


#' Validate transition args
#'
#' Validates the transition args for the samc() function
#'
#' @param x A list
#' @noRd
.validate_tr_args <- function(x) {
  args <- c("fun", "dir", "sym")
  names <- names(x)

  missing_args <- args[!(args %in% names)]
  if (length(missing_args) > 0)
    stop(paste("Missing argument in tr_args:", missing_args), call. = FALSE)

  unknown_args <- names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in tr_args:", unknown_args), call. = FALSE)

  dup_args <- names[duplicated(names)]
  if (length(dup_args) > 0)
    stop(paste("Duplicate argument in tr_args:", dup_args), call. = FALSE)

  if (!is.function(x$fun)) {
    stop("`fun`` must be a function.", call. = FALSE)
  } else if (!(x$dir %in% c(4,8))) {
    stop("`dir` must be set to either 4 or 8", call. = FALSE)
  } else if (!is.logical(x$sym)) {
    stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)
  }
}


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
  signature(samc = "samc", x = "samc_raster"),
  function(samc, x) {

    x = .rasterize(x)

    check(samc, x)

    if (terra::nlyr(x) == 0) {
     stop("Missing absorption data", call. = FALSE)
    }


    abs_mat <- terra::values(x)
    abs_vec <- as.vector(abs_mat[, 1])

    abs_minmax = terra::minmax(x)

    if (min(abs_minmax["min", ]) < 0 || max(abs_minmax["max", ]) > 1) {
      stop("Absorption values must be in range of 0-1", call. = FALSE)
    }

    excl <- which(is.na(abs_vec))
    if (length(excl) > 0) {
      abs_mat <- abs_mat[-excl, , drop = FALSE]
    }

    if (is.null(names(x))) colnames(abs_mat) <- 1:ncol(abs_mat)

    if ("" %in% names(x)) stop("Mix of named and unnamed maps/layers", call. = FALSE)
    if (any(duplicated(names(x)))) stop("Duplicate names", call. = FALSE)

    colnames(abs_mat) <- names(x)

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
  signature(samc = "samc", x = "list"),
  function(samc, x) {

    if (!all(sapply(x, is.matrix))) stop("List can only contain matrices. If using rasters, use raster::stack() or c() to combine terra SpatRasters instead.", call. = FALSE)

    x <- lapply(x, .rasterize)

    if(is.null(names(x))) names(x) = 1:length(x)

    return(.process_abs_states(samc, terra::rast(x)))
  })


#' Process occupancy input
#'
#' Process occupancy input
#'
#' @param samc A samc-class object
#' @param x occupancy input
#' @noRd
setGeneric(
  ".process_occ",
  function(samc, x) {
    standardGeneric(".process_occ")
  })

# TODO: find a way to check the input type for `occ` to the input type to samc()

#' @noRd
setMethod(
  ".process_occ",
  signature(samc = "samc", x = "numeric"),
  function(samc, x) {
    if (any(!is.finite(x)) || any(x < 0)) stop("`occ` input must only contain positive numeric values")

    if (length(x) != nrow(samc$q_matrix)) stop("`occ` input length does not match number of transient states")

    return(x)
  })

#' @noRd
setMethod(
  ".process_occ",
  signature(samc = "samc", x = "SpatRaster"),
  function(samc, x) {
    check(samc, x)

    pv <- as.vector(terra::values(x))
    pv <- pv[is.finite(pv)]

    return(.process_occ(samc, pv))
  })

#' @noRd
setMethod(
  ".process_occ",
  signature(samc = "samc", x = "Raster"),
  function(samc, x) {
    return(.process_occ(samc, .rasterize(x)))
  })

setMethod(
  ".process_occ",
  signature(samc = "samc", x = "matrix"),
  function(samc, x) {
    return(.process_occ(samc, .rasterize(x)))
  })
