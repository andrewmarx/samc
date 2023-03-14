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
  tr = .transition(x, fun, dir, sym)

  nrows = terra::nrow(x)
  ncols = terra::ncol(x)
  cell_nums = terra::cells(x)
  ncells = length(cell_nums)

  edge_counts = Matrix::rowSums(tr > 0)
  edge_counts = tr@p[-1] - tr@p[-length(tr@p)] - 1


  # Stuff for indexing things while looping by cols
  row_offsets = numeric(ncells)
  row_offset_sum = 0

  for (i in 1:ncells) {
    row_offsets[i] = row_offset_sum
    row_offset_sum = row_offset_sum + edge_counts[i]
  }
  row_offsets = row_offsets + 1
  rm(row_offset_sum)

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

  #fidelity = terra::values(fidelity)



  tr_i = tr@i
  tr_p = tr@p
  tr_x = tr@x

  # Fill out CRW lookup
  crw_map = matrix(0L, nrow = sum(edge_counts), ncol = 2)
  row_accesses = numeric(ncells)
  dir_lookup = matrix(c(1:4, NA, 5:8), nrow = 3, byrow = TRUE)

  for (p2 in 1:ncells) {
    #print("")

    p1s = tr_i[(tr_p[p2] + 1) : tr_p[p2 + 1]] + 1
    #print(p1s)

    for (p1 in p1s) {
      if (p1 != p2) { # if an edge
        # crw_map
        src = cell_nums[p1]
        dst = cell_nums[p2]
        dr = dir_lookup[((dst - 1) %/% ncols - (src - 1) %/% ncols) + 2,
                        ((dst - 1) %% ncols - (src - 1) %% ncols) + 2]
        crw_map[row_offsets[p1] + row_accesses[p1], ] =
          c(p1, dr)

        row_accesses[p1] = row_accesses[p1] + 1

      }
    }
  }

  # e2 then e1 for CSC for entire mat
  mat_p = integer(sum(edge_counts) + 1)
  mat_x = numeric(sum((edge_counts + 1) * edge_counts))
  mat_i = integer(sum((edge_counts + 1) * edge_counts))

  sum = 0
  index = 1
  for (p2 in 1:ncells) {
    p3s = tr_i[(tr_p[p2] + 1) : tr_p[p2 + 1]] + 1
    p3i = 0
    for (p3 in p3s) {
      if (p2 != p3) {
        e2i = row_offsets[p2] + p3i

        sum = sum + edge_counts[p2] + 1
        mat_p[e2i + 1] = as.integer(sum)


        p1s = tr_i[(tr_p[p2] + 1) : tr_p[p2 + 1]] + 1
        p1i = 0

        fid_check = FALSE
        fid_offset = 0

        count = 0

        for (p1 in p1s) {
          if (p1 != p2) {

            p2s = tr_i[(tr_p[p1] + 1) : tr_p[p1 + 1]] + 1
            p2s = p2s[p2s != p1]
            p2i = which(p2s == p2)

            e1i = row_offsets[p1] + p2i - 1

            mat_x[index] = tr[p2, p3] * ang_mat[crw_map[e1i, 2], crw_map[e2i, 2]] # TODO the sparse matrix access here is why the function is slow
            mat_i[index] = e1i
          } else {
            mat_x[index] = 0 #fidelity[cell_nums[p1]]
            mat_i[index] = e2i
          }
          index = index + 1
          p1i = p1i + 1
        }

        p3i = p3i + 1

      }
    }

  }
  # mat_p
  # mat_x
  # mat_i

  mat = new("dgCMatrix")
  mat@Dim = c(as.integer(sum(edge_counts)), as.integer(sum(edge_counts)))

  mat@p = mat_p
  mat@i = as.integer(mat_i - 1)
  mat@x = mat_x

  #View(as.matrix(mat))


  # normalization

  fidelity = terra::values(fidelity)
  tmp = 1 - terra::values(absorption) - fidelity
  rs = Matrix::rowSums(mat)

  # TODO figure out where bug is in this optimization
  # crw_lookup = as.vector(crw_map[,1])
  #
  # i_index = 1
  # for (p in 1:sum(edge_counts)) {
  #   row_count = mat_p[p+1] - mat_p[p]
  #   for (i in 1:row_count) {
  #     row = mat_i[i_index] + 1
  #     if (p != row) {
  #       #mat_x[i_index] = i_index # useful for validation
  #       #mat_x[i_index] = cell_nums[row] # useful for validation
  #       #print(c(p, row))
  #       #assign("ts", list(mat_p, mat_i), envir = globalenv())
  #
  #       mat_x[i_index] = -mat_x[i_index]/rs[row] * tmp[cell_nums[crw_lookup[row]]]
  #     } else {
  #       mat_x[i_index] =  1 - fidelity[cell_nums[crw_lookup[row]]]
  #     }
  #     i_index = i_index + 1
  #   }
  # }
  #
  # # For perf reasons
  # mat@x = mat_x


  i_index = 1
  for (p in 1:sum(edge_counts)) {
    row_count = mat@p[p+1] - mat@p[p]
    for (i in 1:row_count) {
      row = mat@i[i_index] + 1
      if (p != row) {
        #mat_x[i_index] = i_index # useful for validation
        #mat_x[i_index] = cell_nums[row] # useful for validation
        #print(c(p, row))
        #assign("ts", list(mat_p, mat_i), envir = globalenv())

        mat@x[i_index] = -mat@x[i_index]/rs[row] * tmp[cell_nums[crw_map[,1][row]]]
      } else {
        mat@x[i_index] =  1 - fidelity[cell_nums[crw_map[,1][row]]]
      }
      i_index = i_index + 1
    }
  }

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
    print(vec)
    print(x)
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
