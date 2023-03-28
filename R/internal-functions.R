# Copyright (c) 2020-2023 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.
#
# This file is for internal functions. They are subject to change and should not
# be used by users.
#

#' RW function
#'
#' An internal function for creating RW samc matrix
#'
#' @noRd
.rw <- function(x, absorption, fidelity, fun, dir, sym) {
  tr = .transition(x, fun, dir, sym)

  rs = Matrix::rowSums(tr)

  tr_i = tr@i
  tr_p = tr@p
  tr_x = tr@x

  tmp = 1 - terra::values(absorption) - terra::values(fidelity)


  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  fidelity = terra::values(fidelity)

  i_index = 1
  for (p in 1:ncells) {
    row_count = tr_p[p+1] - tr_p[p]
    for (i in 1:row_count) {
      row = tr_i[i_index] + 1
      if (p != row) {
        #mat_x[i_index] = i_index # useful for validation
        #mat_x[i_index] = cell_nums[row] # useful for validation
        #print(c(p, row))
        #assign("ts", list(mat_p, mat_i), envir = globalenv())
        tr_x[i_index] = -tr_x[i_index]/rs[row] * tmp[cell_nums[row]]
      } else {
        tr_x[i_index] = 1 - fidelity[cell_nums[row]]
      }
      i_index = i_index + 1
    }
  }

  tr@x = tr_x

  tr
}

#' CRW function
#'
#' An internal function for creating RW samc objects
#'
#' @noRd
.crw <- function(x, absorption, fidelity, fun, dir, sym = TRUE, model) {
  tr = .tr_vals_res(x, dir)
  edge_counts = sum(is.finite(tr))
  edge_nums = tr
  edge_nums[is.finite(edge_nums)] = 1:edge_counts


  nrows = terra::nrow(x)
  ncols = terra::ncol(x)
  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  edge_counts = sum(is.finite(tr))

  # Angle matrix
  # TODO make sure works for 4 directions
  dir_vec = matrix(c(-1, 1,
                      0, 1,
                      1, 1,
                     -1, 0,
                      1, 0,
                     -1, -1,
                      0, -1,
                      1, -1),
                   nrow = 8, byrow = TRUE)

  ang_mat = matrix(nrow = 8, ncol = 8)

  for (r in 1:8) {
    for (c in 1:8) {
      mag_v1 = sqrt(sum(dir_vec[r, ]^2))
      mag_v2 = sqrt(sum(dir_vec[c, ]^2))

      ang_mat[r, c] = circular::dvonmises(circular::circular(acos(sum(dir_vec[r, ] * dir_vec[c, ]) / (mag_v1 * mag_v2))), mu = circular::circular(0), kappa = model$dist$kappa)
    }
  }


  # tr offset vars
  if (dir == 4) {
    dir_vec = c(1:2, 0, 3:4)
    offsets = c(-ncols, -1, 1, ncols) * 4
  } else if (dir == 8) {
    dir_vec = c(1:4, 0, 5:8)
    offsets = c(-ncols - 1, -ncols, -ncols + 1, -1, 1, ncols - 1, ncols, ncols + 1) * 8
  }


  # Fill out mat_p
  mat_p = integer(edge_counts + 1)
  cell = 0

  for (r in 1:nrows) {
    vals = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
    for (c in 1:ncols) {
      cell = cell + 1
      if (is.finite(vals[c])) {
        p1 = (cell - 1) * dir

        # loop through valid edges
        for (d in 1:dir) {
          e1 = p1 + d

          if (!is.na(tr[e1])) {
            p2 = p1 + offsets[d]

            # loop through valid secondary edges
            for (dv in dir_vec) {
              if (dv) {
                e2 = p2 + dv
                if (!is.na(tr[e2])) {
                  mat_p[edge_nums[e2] + 1] = mat_p[edge_nums[e2] + 1] + 1
                }
              } else {
                mat_p[edge_nums[e1] + 1] = mat_p[edge_nums[e1] + 1] + 1
              }
            }
          }
        }
      }
    }
  }

  for (i in 2:length(mat_p)) {
    mat_p[i] = mat_p[i] + mat_p[i - 1]
  }



  mat_p_count = integer(edge_counts)

  mat_x = numeric(mat_p[length(mat_p)])
  mat_i = integer(mat_p[length(mat_p)])

  crw_map = matrix(0L, nrow = edge_counts, ncol = 2)


  cell = 0
  index = 0

  # Loop through cells with values
  for (r in 1:nrows) {
    fid = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
    vals = 1 - fid - terra::values(absorption, mat = FALSE, row = r, nrows = 1)
    for (c in 1:ncols) {
      cell = cell + 1
      if (is.finite(vals[c])) {
        p1 = (cell - 1) * dir

        # loop through valid edges
        for (d in 1:dir) {
          e1 = p1 + d
          if (!is.na(tr[e1])) {
            p2 = p1 + offsets[d]
            e1_num = edge_nums[e1]

            index = index + 1
            crw_map[index, ] = c(cell, d)

            rs = 0
            fid_index = NA

            row_indices = numeric(dir)
            row_indices[] = NA

            # loop through valid secondary edges
            for (dv in dir_vec) {
              if (dv) {
                e2 = p2 + dv

                if (!is.na(tr[e2])) {
                  e2_num = edge_nums[e2]
                  mat_p_count[e2_num] = mat_p_count[e2_num] + 1

                  res = tr[e2] * ang_mat[d, dv]
                  rs = rs + res

                  row_indices[dv] = mat_p[e2_num] + mat_p_count[e2_num]
                  mat_x[mat_p[e2_num] + mat_p_count[e2_num]] = -res * vals[c]
                  mat_i[mat_p[e2_num] + mat_p_count[e2_num]] = e1_num
                }
              } else {
                mat_p_count[e1_num] = mat_p_count[e1_num] + 1

                fid_index = mat_p[e1_num] + mat_p_count[e1_num]
                mat_x[fid_index] = 1 - fid[c]
                mat_i[fid_index] = e1_num
              }
            }
            row_indices = row_indices[!is.na(row_indices)]
            mat_x[row_indices] = mat_x[row_indices]/rs
          }
        }
      }
    }
  }

  mat_i = mat_i - 1

  mat = new("dgCMatrix")
  mat@Dim = c(as.integer(sum(edge_counts)), as.integer(sum(edge_counts)))

  mat@p = as.integer(mat_p)
  mat@i = as.integer(mat_i)
  mat@x = mat_x

  #View(as.matrix(mat))

  return(
    list(tr = mat,
         crw = crw_map,
         abs = terra::values(absorption)[cell_nums[crw_map[,1]]])
  )
}

#' Transition function
#'
#' An internal function for creating transition matrices
#'
#' @noRd
.transition <- function(x, fun, dir, sym = TRUE) {

  lonlat = terra::is.lonlat(x)

  #data_crs = terra::crs(resistance)
  #ncells = terra::ncell(resistance)
  nrows = terra::nrow(x)
  ncols = terra::ncol(x)

  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  cell_lookup = matrix(0, ncols, nrows) # Intentionally reversed for rast->mat
  cell_lookup[cell_nums] = 1:length(cell_nums)

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

  n_pair = n_pair*2 + ncells

  mat = new("dgCMatrix")
  mat@Dim = c(ncells, ncells)



  if (dir == 4) {
    offset_r = c(0, -1, 1, 0)
    offset_c = c(-1, 0, 0, 1)

    offset_n = c(-ncols, -1, 1, ncols)
  } else if (dir == 8) {
    offset_r = c(-1, 0, 1, -1, 1, -1, 0, 1)
    offset_c = c(-1, -1, -1, 0, 0, 1, 1, 1)

    offset_n = c(-ncols-1, -ncols, -ncols+1, -1, 1, ncols - 1, ncols, ncols+1) # ncols of raster = nrow() of matrix
  } else {
    stop("Issue with directions", call. = FALSE)
  }

  dir1 = 1:(dir/2)
  dir2 = dir1 + (dir/2)

  mat_p = integer(ncells+1)

  mat_i = integer(n_pair)
  i_index = 1

  mat_x = numeric(n_pair)

  row_count = 0L


  dist = .build_lookup_function(x, dir)

  nc = nrows
  nr = ncols

  x = terra::values(x)

  for (i in 1:length(cell_nums)) {
    num = cell_nums[i]

    row = (num - 1) %% nr + 1
    col = (num - 1) %/% nr + 1

    #print(paste(row, col))
    rows = row + offset_r
    cols = col + offset_c
    nums = num + offset_n

    rc = !(rows < 1 | rows > nr | cols < 1 | cols > nc)

    for (d in dir1) {
      #print(d)
      if (rc[d]) {
        mat_row = cell_lookup[nums[d]]

        if (mat_row) {
          row_count = row_count + 1L

          result = fun(x[c(num, nums[d])]) / dist(num, d)

          mat_i[i_index] = as.integer(mat_row - 1)
          mat_x[i_index] = result

          i_index = i_index + 1
        }
      }
    }


    row_count = row_count + 1L

    mat_row = cell_lookup[num]

    mat_i[i_index] = as.integer(mat_row - 1)
    mat_x[i_index] = 0

    i_index = i_index + 1


    for (d in dir2) {
      if (rc[d]) {
        mat_row = cell_lookup[nums[d]]

        if (mat_row) {
          row_count = row_count + 1L

          result = fun(x[c(num, nums[d])]) / dist(num, d)

          mat_i[i_index] = as.integer(mat_row - 1)
          mat_x[i_index] = result

          i_index = i_index + 1
        }
      }
    }

    mat_p[i+1] = row_count
  }

  mat@p = mat_p
  mat@i = mat_i
  mat@x = mat_x

  return(mat)
}


.tr_vals = function(data, fun, dir) {

  nrows = terra::nrow(data)
  ncols = terra::ncol(data)

  result = numeric(nrows * ncols * dir)
  index = 0

  dist = .build_lookup_mat(data, dir)

  if (dir == 4) {
    dir = c(2, 4, 6, 8)
  } else if (dir == 8) {
    dir = c(1:4, 6:9)
  }

  for (r in 1:nrows) {
    vals = terra::focalValues(data, 3, r,1)

    for (c in 1:ncols) {
      v = vals[c, 5]
      for (d in dir) {
        index = index + 1

        result[index] = fun(c(v, vals[c, d])) / dist[r, d]

      }
    }
  }

  result
}


.tr_vals_res = function(data, dir) {

  nrows = terra::nrow(data)
  ncols = terra::ncol(data)

  result = numeric(nrows * ncols * dir)
  index = 0

  dist = .build_lookup_mat(data, dir)

  if (dir == 4) {
    dir = c(2, 4, 6, 8)
  } else if (dir == 8) {
    dir = c(1:4, 6:9)
  }

  for (r in 1:nrows) {
    vals = terra::focalValues(data, 3, r,1)

    for (c in 1:ncols) {
      v = vals[c, 5]
      for (d in dir) {
        index = index + 1

        result[index] = 2 / ((v + vals[c, d]) * dist[r, d])

      }
    }
  }

  result
}


#' Build direction lookup function
#'
#' TODO description here
#'
#' @param x A function
#' @noRd

.build_lookup_function <- function(rast, dir) {

  # TODO Create tests to directly validate results for both directions options for planar and latlon

  lonlat = terra::is.lonlat(rast)
  nrows = terra::nrow(rast)
  ncols = terra::ncol(rast)

  if (dir == 4) {
    dist_lookup = c(1, 1, 1, 1)
  } else if (dir == 8) {
    dist_lookup = c(sqrt(2), 1, sqrt(2), 1, 1, sqrt(2), 1, sqrt(2))
  }

  dist = function(x, dir) {
    # x not used intentionally
    dist_lookup[dir]
  }

  if(!is.na(lonlat)) {
    if (lonlat) {
      cn = (0:(nrows - 1)) * ncols + 1

      adj = terra::adjacent(rast, cn, directions = 8, pairs = TRUE)

      dist = terra::distance(
        terra::xyFromCell(rast, adj[, 1]),
        terra::xyFromCell(rast, adj[, 2]),
        lonlat = TRUE, pairwise = TRUE)

      adj = terra::adjacent(rast, cn, directions = dir, pairs = FALSE)
      adj = t(adj)

      dist_lookup = adj
      dist_lookup[!is.nan(dist_lookup)] = dist

      if (dir == 4) {
        dir_reindex = c(1, 3, 3, 4)
      } else if (dir == 8) {
        dir_reindex = c(3, 2, 3, 5, 5, 8, 7, 8)
      } else {
        stop("Invalid directions", call. = FALSE)
      }


      dist <- function(x, dir) {
        dir = dir_reindex[dir]

        x = trunc((x - 1) / ncols) + 1

        dist_lookup[dir, x]
      }
    }
  }

  return(dist)
}

.build_lookup_mat <- function(rast, dir) {

  # TODO Create tests to directly validate results for both directions options for planar and latlon

  lonlat = terra::is.lonlat(rast)
  nrows = terra::nrow(rast)
  ncols = terra::ncol(rast)

  mat = matrix(nrow = nrows, ncol = 9)
  mat[, c(2, 4, 6, 8)] = 1

  if (dir == 8) {
    mat[, c(1, 3, 7, 9)] = sqrt(2)
  }


  if(!is.na(lonlat)) {
    if (lonlat) {
      cn = (0:(nrows - 1)) * ncols + 1

      adj = terra::adjacent(rast, cn, directions = 8, pairs = TRUE)

      dist = terra::distance(
        terra::xyFromCell(rast, adj[, 1]),
        terra::xyFromCell(rast, adj[, 2]),
        lonlat = TRUE, pairwise = TRUE)

      adj = terra::adjacent(rast, cn, directions = dir, pairs = FALSE)
      adj = t(adj)

      adj[!is.nan(adj)] = dist
      adj = t(adj)

      mat[, c(1:4, 6:9)] = adj
      mat[, 1] = mat[, 3]
      mat[, 4] = mat[, 6]
      mat[, 7] = mat[, 9]
    }
  }

  return(mat)
}

#' Build turning lookup function
#'
#' TODO description here
#'
#' @param x A function
#' @noRd

.build_turn_function <- function() {

}


#' Convolution
#'
#' TODO description here
#'
#' @param x A function
#' @noRd

.convolution <- function(res, abso, fid, dir, sym, threads) {

  if (dir == 4) {
    kernel = matrix(
      c(0, 1, 0,
        1, 0, 1,
        0, 1, 0), 3)
  } else if (dir == 8) {
    kernel = matrix(
      c(1 / sqrt(2), 1, 1 / sqrt(2),
        1.0000000, 0, 1.0000000,
        1 / sqrt(2), 1, 1 / sqrt(2)
      ), 3, 3
    )
  } else {
    stop("'directions' must be equal to either 4 or 8")
  }

  nr = nrow(res)
  nc = ncol(res)

  res = as.matrix(res)
  abso = as.matrix(abso)
  fid = as.matrix(fid)

  res[!is.finite(res)] = 0
  abso[!is.finite(abso)] = 0
  fid[!is.finite(fid)] = 0

  dim(res) = c(nc, nr)
  dim(abso) = c(nc, nr)
  dim(fid) = c(nc, nr)

  samc:::.build_convolution_cache(kernel, res, fid, abso, sym, threads)
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

  # if (any(x > 10000))
  #   stop("Due to how the short-term metrics are calculated and the way that
  # decimal numbers are handled by computers, numerical issues related to
  # precision arise when a time step value is too high. Currently, a hard
  # limit of 10000 time steps is enforced to encourage users to more
  # seriously consider how many time steps are relevant to their use case.
  # For example, if a single time step represents 1 day, then the current
  # limit represents 24.7 years. There is flexibility to increase the limit
  # if a justification can be made for it, but it's far more likely that
  # users will generally want far fewer time steps for ecologically relevant
  # results and to avoid the cummulative precision issues.", call. = FALSE)
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
    #print(vec)
    #print(x)
    stop(paste("\nInvalid location name:", invalid_names), call. = FALSE)
  }
}


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

    if (samc@model$name == "CRW") {
      stop("CRW model requires a list with location and direction.", call. = FALSE)
    } else {
      .validate_locations(samc, x)
    }
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
.validate_model <- function(x, method) {
  args <- c("name", "fun", "dir", "sym")
  names <- names(x)

  dup_args <- names[duplicated(names)]
  if (length(dup_args) > 0)
    stop(paste("Duplicate argument in model:", dup_args), call. = FALSE)

  if (!("name" %in% names)) {
    x$name = "RW"
    names = c(names, "name")
  }

  if (method == "conv") {
    if (x$name != "RW") stop("Convolution only supports RW model.")
  }


  if (x$name == "CRW") {
    args = c(args, "dist")
  }

  missing_args <- args[!(args %in% names)]
  if (length(missing_args) > 0)
    stop(paste("Missing argument in model:", missing_args), call. = FALSE)

  unknown_args <- names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)

  if (!is.function(x$fun)) {
    stop("`fun`` must be a function.", call. = FALSE)
  } else if (!(x$dir %in% c(4,8))) {
    stop("`dir` must be set to either 4 or 8", call. = FALSE)
  } else if (!is.logical(x$sym)) {
    stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)
  }

  if (x$name == "CRW") {
    dist = x$dist

    dist_args = c("name")
    dist_names = names(dist)

    if (!("name" %in% dist_names)) {
      stop("Missing distribution name.", call. = FALSE)
    } else if (dist$name == "vonMises") {
      dist_args = c(dist_args, "kappa")
    } else {
      stop(paste("Invalid distribution name:", dist$name), call. = FALSE)
    }

    missing_args <- dist_args[!(dist_args %in% dist_names)]
    if (length(missing_args) > 0)
      stop(paste("Missing argument in dist:", missing_args), call. = FALSE)

    unknown_args <- dist_names[!(dist_names %in% dist_args)]
    if (length(unknown_args) > 0)
      stop(paste("Unknown argument in dist:", unknown_args), call. = FALSE)

    if (dist$name == "vonMises") {
      if (!is.numeric(dist$kappa))
        stop("kappa must be single non-negative numeric value.", call. = FALSE)

      if (length(dist$kappa) != 1)
        stop("kappa must be single non-negative numeric value.", call. = FALSE)

      if (!is.finite(dist$kappa))
        stop("kappa must be single non-negative numeric value.", call. = FALSE)

      if (dist$kappa < 0)
        stop("kappa must be single non-negative numeric value.", call. = FALSE)
    }
  }

  return(x)
}


#' Validate options
#'
#' Validates the options args for the samc() function
#'
#' @param x A list
#' @noRd
.validate_options <- function(x) {
  if (is.null(x)) {
    x = list(
      method = "direct",
      threads = 1,
      override = FALSE
    )
  } else if (is.list(x)) {
    # TODO finish this
  } else {
    stop("options argument must be a list or left empty for default values", call. = FALSE)
  }

  x
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
  signature(samc = "samc", x = "SpatRaster"),
  function(samc, x) {

    x = rasterize(x)

    check(samc@map, x)

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

    x <- lapply(x, rasterize)

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

    if (length(x) != length(samc@data@t_abs)) stop("`init` input length does not match number of transient states")

    return(x)
  })

#' @noRd
setMethod(
  ".process_init",
  signature(samc = "samc", x = "SpatRaster"),
  function(samc, x) {
    check(samc@map, x)

    pv <- as.vector(terra::values(x))
    pv <- pv[is.finite(pv)]

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
