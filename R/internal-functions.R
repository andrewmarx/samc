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
