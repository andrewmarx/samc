# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.


#' RW function
#'
#' An internal function for creating RW samc matrix
#'
#' @noRd
.rw <- function(x, absorption, fidelity, model) {
  dir = model$dir
  fun = model$fun
  sym = model$sym

  if (model$name == "RW") {
    tr = .tr_vals(x, fun, dir)
  } else if (model$name == "SSF") {
    tr = .tr_vals_ssf(x, fun, dir, model$ssc)
  } else {
    stop("", call. = FALSE)
  }

  nedges = sum(is.finite(tr))

  nrows = terra::nrow(x)
  ncols = terra::ncol(x)
  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  cell_lookup = matrix(0, ncols, nrows) # Intentionally reversed for rast->mat
  cell_lookup[cell_nums] = 1:length(cell_nums)

  # tr offset vars
  if (dir == 4) {
    dir_vec = c(1:2, 0, 3:4)
    offsets = c(-ncols, -1, 1, ncols)
  } else if (dir == 8) {
    dir_vec = c(1:4, 0, 5:8)
    offsets = c(-ncols - 1, -ncols, -ncols + 1, -1, 1, ncols - 1, ncols, ncols + 1)
  }

  # Fill out mat_p
  mat_p = integer(ncells + 1)
  mat_p[] = 1 # every column will have a fidelity value
  mat_p[1] = 0 # except first entry of mat_p does not refer to a column
  cell = 0

  for (r in 1:nrows) {
    vals = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
    for (c in 1:ncols) {
      cell = cell + 1
      if (is.finite(vals[c])) {
        p1 = cell_lookup[cell]
        p1i = (cell - 1) * dir

        # loop through valid edges
        for (d in dir_vec) {
          if (d) {
            p2i = p1i + d

            if (!is.na(tr[p2i])) {
              p2 = cell_lookup[cell + offsets[d]]

              mat_p[p2 + 1] = mat_p[p2 + 1] + 1
            }
          }
        }
      }
    }
  }

  for (i in 2:length(mat_p)) {
    mat_p[i] = mat_p[i] + mat_p[i - 1]
  }

  mat_p_count = integer(ncells + 1)

  mat_x = numeric(nedges + ncells)
  mat_i = integer(nedges + ncells)

  cell = 0

  row_indices = numeric(dir)

  for (r in 1:nrows) {
    fid = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
    vals = 1 - fid - terra::values(absorption, mat = FALSE, row = r, nrows = 1)

    for (c in 1:ncols) {
      cell = cell + 1
      if (is.finite(vals[c])) {
        p1 = cell_lookup[cell]
        p1i = (cell - 1) * dir

        # loop through valid edges
        rs = 0
        fid_index = NA

        row_indices[] = NA

        for (d in dir_vec) {
          if (d) {
            p2i = p1i + d

            if (!is.na(tr[p2i])) {
              p2 = cell_lookup[cell + offsets[d]]

              mat_p_count[p2] = mat_p_count[p2] + 1

              res = tr[p2i]
              rs = rs + res

              row_indices[d] = mat_p[p2] + mat_p_count[p2]
              mat_x[mat_p[p2] + mat_p_count[p2]] = -res * vals[c]
              mat_i[mat_p[p2] + mat_p_count[p2]] = p1
            }
          } else {
            mat_p_count[p1] = mat_p_count[p1] + 1

            fid_index = mat_p[p1] + mat_p_count[p1]
            mat_x[fid_index] = 1 - fid[c]
            mat_i[fid_index] = p1
          }
        }
        row_indices = row_indices[!is.na(row_indices)]
        mat_x[row_indices] = mat_x[row_indices]/rs
      }
    }
  }

  mat_i = mat_i - 1

  mat = new("dgCMatrix")
  mat@Dim = c(as.integer(ncells), as.integer(ncells))

  mat@p = as.integer(mat_p)
  mat@i = as.integer(mat_i)
  mat@x = mat_x

  mat
}


#' CRW function
#'
#' An internal function for creating RW samc objects
#'
#' @noRd
.crw <- function(x, absorption, fidelity, fun, dir, sym = TRUE, model) {
  tr = .tr_vals(x, fun, dir)

  edge_counts = sum(is.finite(tr))
  edge_nums = tr
  edge_nums[is.finite(edge_nums)] = 1:edge_counts


  nrows = terra::nrow(x)
  ncols = terra::ncol(x)
  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  edge_counts = sum(is.finite(tr))

  mu = circular::circular(0)

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

      ang_mat[r, c] = circular::circular(acos(sum(dir_vec[r, ] * dir_vec[c, ]) / (mag_v1 * mag_v2)))
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
  crw_index = 0
  index = 0

  # Loop through cells with values
  for (r in 1:nrows) {
    fid = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
    vals = 1 - fid - terra::values(absorption, mat = FALSE, row = r, nrows = 1)
    kappa_vals = terra::values(model$kappa, mat = FALSE, row = r, nrows = 1)
    for (c in 1:ncols) {
      cell = cell + 1
      if (is.finite(vals[c])) {
        p1 = (cell - 1) * dir
        index = index + 1

        # loop through valid edges
        for (d in 1:dir) {
          e1 = p1 + d
          if (!is.na(tr[e1])) {
            p2 = p1 + offsets[d]
            e1_num = edge_nums[e1]

            crw_index = crw_index + 1
            crw_map[crw_index, ] = c(index, d)

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

                  res = tr[e2] * circular::dvonmises(ang_mat[d, dv], mu = mu, kappa = kappa_vals[c])
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


  dim(tr) = c(dir, length(tr)/dir)
  tr = scale(tr, FALSE, colSums(tr, na.rm = TRUE))
  attr(tr, 'scaled:scale') = NULL

  return(
    list(tr = mat,
         crw = crw_map,
         prob = tr,
         abs = terra::values(absorption)[cell_nums[crw_map[,1]]])
  )
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

  if (is(fun, "function")) {
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
  } else if (fun == "1/mean(x)") {
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
  } else {
    stop("Invalid transition function defined", call. = FALSE)
  }

  result
}


.tr_vals_ssf = function(data, fun, dir, ssc) {

  nrows = terra::nrow(data)
  ncols = terra::ncol(data)

  result = numeric(nrows * ncols * dir)
  index = 0

  dist = .build_lookup_mat(data, dir)^(ssc)

  if (dir == 4) {
    dir = c(2, 4, 6, 8)
  } else if (dir == 8) {
    dir = c(1:4, 6:9)
  }

  if (is(fun, "function")) {
    for (r in 1:nrows) {
      vals = terra::focalValues(data, 3, r,1)

      for (c in 1:ncols) {
        v = vals[c, 5]
        for (d in dir) {
          index = index + 1

          result[index] = fun(c(v, vals[c, d])) * dist[r, d]
        }
      }
    }
  } else if (fun == "x[2]") {
    for (r in 1:nrows) {
      vals = terra::focalValues(data, 3, r,1)

      for (c in 1:ncols) {
        for (d in dir) {
          index = index + 1

          result[index] = vals[c, d] * dist[r, d]
        }
      }
    }
  } else {
    stop("Invalid transition function defined", call. = FALSE)
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

  .build_convolution_cache(kernel, res, fid, abso, sym, threads)
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
  # results and to avoid the cumulative precision issues.", call. = FALSE)
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

  if (any(x > samc@nodes))
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


#' Map Location
#'
#' Map a location
#'
#' @param x A list
#' @noRd
.map_location <- function(samc, x) {
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

#' Used to disable CRW
#'
#' Disable CRW
#'
#' @param samc samc model
#' @noRd
.disable_crw <- function(samc) {
  if (samc@model$name == "CRW") stop("Metric/parameter combination not currently supported for CRW", call. = FALSE)
}

#' Used to disable convolution
#'
#' Disable convolution
#'
#' @param samc samc model
#' @noRd
.disable_conv <- function(samc) {
  if (samc@solver == "conv") stop("Metric/parameter combinaton not currently supported for the convolution algorithm", call. = FALSE)
}


#' Summarize CRW
#'
#' Summarize CRW
#'
#' @param samc samc model
#' @noRd
.summarize_crw <- function(samc, vec, fun) {
  if (length(vec) != nrow(samc@crw_map))
    stop("The length of the vector does not match the number of non-NA cells in the landscape data", call. = FALSE)

  aggregate(vec ~ samc@crw_map[, 1], FUN = fun)$vec
}
