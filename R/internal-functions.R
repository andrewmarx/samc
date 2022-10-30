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

  fid_max = terra::minmax(fidelity)[2]
  if (fid_max == 0) {
    fidelity = 0
  }

  cell_nums = terra::cells(resistance)
  adj = terra::adjacent(resistance, cells=cell_nums, pairs=TRUE, directions=dir)
  adj = adj[adj[, 2] %in% cell_nums, ]

  adj_length = nrow(adj)

  if(sym) adj = adj[adj[, 1] < adj[, 2], ]

  coords <- cbind(terra::xyFromCell(resistance, adj[, 1]),
                  terra::xyFromCell(resistance, adj[, 2]))

  i = adj[,1]
  j = adj[,2]

  sums = numeric(max(adj))

  resistance = terra::values(resistance)
  dim(resistance) = nrow(resistance)

  adj[] = resistance[adj]
  rm(resistance);gc()

  dist <- terra::distance(coords[,1:2],
                          coords[,3:4],
                          lonlat = terra::is.lonlat(absorption),
                          pairwise = TRUE)

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

  transition.values = transition.values / sums[i]
  rm(sums); gc()

  tmp = as.vector(terra::values(1 - absorption - fidelity))
  transition.values = tmp[i] * transition.values
  rm(tmp); gc()

  if(!all(transition.values >= 0)){
    warning("transition function gives negative values")
  }

  if (fid_max > 0) {
    i = c(i, cell_nums)
    j = c(j, cell_nums)
    transition.values = c(transition.values, as.vector(terra::values(fidelity))[cell_nums])
  }


  # Adjust for NAs
  i = match(i, cell_nums)
  j = match(j, cell_nums)


  transitionMatrix = Matrix::sparseMatrix(i = i, j = j, x = transition.values)

  return(transitionMatrix)
}


.transition_new <- function(resistance, absorption, fidelity, fun, dir, sym = TRUE) {
  if (class(fun) == "character" || !(dir %in% c(4, 8))) {
    stop("gdistance's named funtion options not supported")
    #return(gdistance::transition(data, fun, dir, sym))
  }

  lonlat = terra::is.lonlat(absorption)


  #data_crs = terra::crs(resistance)
  #ncells = terra::ncell(resistance)
  nrows = terra::nrow(resistance)
  ncols = terra::ncol(resistance)

  fid_minmax = terra::minmax(fidelity)
  fid_max = fid_minmax["max", 1]
  if (fid_max == 0) {
    fidelity = 0
  }

  cell_nums = terra::cells(resistance)
  ncells = length(cell_nums)

  cell_lookup = numeric(ncells) # Pack down to bits and use bitwise operators
  cell_lookup[cell_nums] = cell_nums
  dim(cell_lookup) = c(ncols, nrows)

  n_pair = 0


  for (i in 1:(nrow(cell_lookup) - 1)) {
    n_pair = n_pair + sum(cell_lookup[i, ] & cell_lookup[i+1, ])
  }
  gc()

  for (i in 1:(ncol(cell_lookup) - 1)) {
    n_pair = n_pair + sum(cell_lookup[, i] & cell_lookup[, i+1])
  }
  gc()

  if (dir == 8) {
    for (i in 1:(ncol(cell_lookup) - 1)) {
      n_pair = n_pair + sum(cell_lookup[-nrow(cell_lookup), i] & cell_lookup[-1, i + 1])
      n_pair = n_pair + sum(cell_lookup[-1, i] & cell_lookup[-nrow(cell_lookup), i + 1])
    }
  }
  gc()

  #rm(cell_lookup);gc()

  n_pair = n_pair*2

  if (fid_max > 0) {
    fidelity = terra::values(fidelity)
    n_pair = n_pair + sum(fidelity[cell_nums] > 0)
  }

  #rm(cell_nums);gc()


  mat = new("dgCMatrix")
  mat@Dim = c(ncells, ncells)

  #mat@x = 1:n_pair
  #mat@i = 1:n_pair
  #mat@p = 1:(ncells+1)


  resistance = terra::values(resistance)
  dim(resistance) = c(ncols, nrows)

  if (dir == 4) {
    offset_r = c(-1, 0, 0, 1)
    offset_c = c(0, -1, 1, 0)

    offset_n = c(-ncols, -1, 1, ncols)
  } else if (dir == 8) {
    offset_r = c(-1, 0, 1, -1, 1, -1, 0, 1)
    offset_c = c(-1, -1, -1, 0, 0, 1, 1, 1)

    #    offset_r = c(-1, -1, -1, 0, 0, 1, 1, 1)
    #    offset_c = c(-1, 0, 1, -1, 1, -1, 0, 1)

    #offset_n = c(-ncols-1, -1, ncols-1, -ncols, ncols, -ncols+1, 1, ncols+1)
    offset_n = c(-ncols-1, -ncols, -ncols+1, -1, 1, ncols - 1, ncols, ncols+1) # ncols of raster = nrow() of matrix
  } else {
    stop("Issue with directions", call. = FALSE)
  }

  dir1 = 1:(dir/2)
  dir2 = dir1 + (dir/2)

  mat_p = integer(ncells+1)

  mat_i = integer(n_pair)
  i_index = 1
  #i_lookup = match(cell_nums)

  mat_x = numeric(n_pair)

  row_sum = numeric(ncells)


  row_count = 0L

  if (lonlat) {
    warning("geocorrection for latlon not implemented", call. = FALSE)
    dist <- function(r1, c1, r2, c2) {
      1 # TODO update
    }
  } else {
    dist_lookup = c(1, sqrt(2))
    dist = function(r1, c1, r2, c2) {
      dist_lookup[abs(r2 - r1) + abs(c2 - c1)]
    }
  }

  for (i in 1:length(cell_nums)) {
    num = cell_nums[i]

    row = (num - 1) %% nrow(resistance) + 1
    col = (num - 1) %/% nrow(resistance) + 1

    #print(paste(row, col))
    rows = row + offset_r
    cols = col + offset_c
    nums = num + offset_n

    #xy = terra::xyFromCell(absorption, num)
    #offset_xy = terra::xyFromCell(absorption, nums)

    #dists = terra::distance(terra::xyFromCell(absorption, num), terra::xyFromCell(absorption, nums), lonlat)

    rc = !(rows < 1 | rows > nrow(resistance) | cols < 1 | cols > ncol(resistance))

    #if (i==2) return(list(num, row, col, rows, cols, rc, nrow(resistance), ncol(resistance)))
    for (d in dir1) {
      #print(d)
      if (rc[d]) {
        row_count = row_count + 1L
        #assign("temp", list(i, d, nums, nums[d], cell_lookup[nums[d]], row, col, rows, cols), globalenv())

        mat_row = cell_lookup[nums[d]]
        result = fun(resistance[c(num, nums[d])]) / dist(row, col, rows[d], cols[d])

        mat_i[i_index] = as.integer(mat_row - 1)
        mat_x[i_index] = result

        row_sum[mat_row] = row_sum[mat_row] + result

        i_index = i_index + 1
      }
    }

    if (fid_max > 0) {
      if (fidelity[num] > 0) {
        row_count = row_count + 1L
        #assign("temp", list(i, d, nums, nums[d], cell_lookup[nums[d]], row, col, rows, cols), globalenv())

        mat_row = cell_lookup[num]

        mat_i[i_index] = as.integer(mat_row - 1)
        mat_x[i_index] = fidelity[num]

        #row_sum[mat_row] = row_sum[mat_row] + fidelity[num] not needed I think

        i_index = i_index + 1
      }
    }

    for (d in dir2) {
      if (rc[d]) {
        #print(d)
        row_count = row_count + 1L
        #assign("temp", list(i, d, nums, nums[d], cell_lookup[nums[d]], row, col, rows, cols), globalenv())

        mat_row = cell_lookup[nums[d]]
        result = fun(resistance[c(num, nums[d])]) / dist(row, col, rows[d], cols[d])

        mat_i[i_index] = as.integer(mat_row - 1)
        mat_x[i_index] = result

        row_sum[mat_row] = row_sum[mat_row] + result

        i_index = i_index + 1
      }
    }
    mat_p[i+1] = row_count
  }

  tmp = terra::values(1 - absorption) - fidelity

  i_index = 1
  for (p in 1:ncells) {
    row_count = mat_p[p+1] - mat_p[p]
    for (i in 1:row_count) {
      row = mat_i[i_index] + 1
      if (p != row) {
        #mat_x[i_index] = i_index # useful for validation
        #mat_x[i_index] = cell_nums[row] # useful for validation
        mat_x[i_index] = mat_x[i_index]/row_sum[row] * tmp[cell_nums[row]]
      }
      i_index = i_index + 1
    }
  }

  mat@p = mat_p
  mat@i = mat_i
  mat@x = mat_x

  return(mat)
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
