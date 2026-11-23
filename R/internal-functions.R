# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.


#' Neighbor cell offsets
#'
#' Linear cell index offsets from a cell to its neighbors, in the same direction
#' order used by `.tr_vals()`. Returns an `ncols` x `dir` matrix so that the
#' offsets can differ by column: on a global lonlat raster, `terra::focalValues()`
#' wraps columns across the antimeridian, so the W neighbors of the first column
#' are in the last column (same row), and vice versa for E neighbors.
#'
#' @noRd
.nb_offsets = function(x, dir) {
    ncols = terra::ncol(x)

    if (dir == 4) {
        dr = c(-1, 0, 0, 1)
        dc = c(0, -1, 1, 0)
    } else if (dir == 8) {
        dr = c(-1, -1, -1, 0, 0, 1, 1, 1)
        dc = c(-1, 0, 1, -1, 1, -1, 0, 1)
    } else {
        stop("Invalid directions", call. = FALSE)
    }

    offsets = matrix(dr * ncols + dc, nrow = ncols, ncol = dir, byrow = TRUE)

    if (isTRUE(terra::is.lonlat(x, global = TRUE))) {
        offsets[1, dc == -1] = offsets[1, dc == -1] + ncols
        offsets[ncols, dc == 1] = offsets[ncols, dc == 1] - ncols
    }

    offsets
}


#' Direction order for matrix entries
#'
#' The order in which a cell's entries are visited when filling the matrix:
#' the direction indices used by `.tr_vals()`, with a 0 marking where the
#' diagonal (fidelity) entry falls so that row indices stay sorted within each
#' column.
#'
#' @noRd
.dir_order = function(dir) {
    if (dir == 4) {
        c(1:2, 0, 3:4)
    } else if (dir == 8) {
        c(1:4, 0, 5:8)
    } else {
        stop("Invalid directions", call. = FALSE)
    }
}


#' Build a dgCMatrix
#'
#' Create an `n` x `n` dgCMatrix directly from its slots
#'
#' @param p Integer vector of column pointers (0-based, length `n + 1`)
#' @param i Integer vector of row indices (0-based)
#' @param x Numeric vector of values
#' @param n Number of rows/columns
#' @noRd
.csc = function(p, i, x, n) {
    mat = new("dgCMatrix")
    mat@Dim = c(as.integer(n), as.integer(n))

    mat@p = p
    mat@i = i
    mat@x = x

    mat
}


#' RW function
#'
#' An internal function for creating RW samc matrix
#'
#' @noRd
.rw = function(x, absorption, fidelity, model) {
    dir = model$dir
    fun = model$fun

    if (model$name == "RW") {
        tr = .tr_vals(x, fun, dir)
    } else if (model$name == "SSF") {
        tr = .tr_vals(x, fun, dir, ssc = model$ssc)
    } else {
        stop("", call. = FALSE)
    }

    nedges = sum(is.finite(tr))

    nrows = terra::nrow(x)
    ncols = terra::ncol(x)
    cell_nums = terra::cells(x)
    ncells = length(cell_nums)

    # Integer so that values read from it don't coerce mat_p/mat_i to double
    cell_lookup = integer(nrows * ncols)
    cell_lookup[cell_nums] = seq_along(cell_nums)

    # tr offset vars
    dir_vec = .dir_order(dir)
    offsets = .nb_offsets(x, dir)

    # Fill out mat_p
    mat_p = integer(ncells + 1)
    mat_p[] = 1L # every column will have a fidelity value
    mat_p[1] = 0L # except first entry of mat_p does not refer to a column
    cell = 0

    for (r in 1:nrows) {
        vals = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
        for (c in 1:ncols) {
            cell = cell + 1
            if (is.finite(vals[c])) {
                p1i = (cell - 1) * dir

                # loop through valid edges
                for (d in dir_vec) {
                    if (d) {
                        p2i = p1i + d

                        if (!is.na(tr[p2i])) {
                            p2 = cell_lookup[cell + offsets[c, d]]

                            mat_p[p2 + 1L] = mat_p[p2 + 1L] + 1L
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

                row_indices[] = NA

                for (d in dir_vec) {
                    if (d) {
                        p2i = p1i + d

                        if (!is.na(tr[p2i])) {
                            p2 = cell_lookup[cell + offsets[c, d]]

                            mat_p_count[p2] = mat_p_count[p2] + 1L

                            res = tr[p2i]
                            rs = rs + res

                            row_indices[d] = mat_p[p2] + mat_p_count[p2]
                            mat_x[mat_p[p2] + mat_p_count[p2]] = -res * vals[c]
                            mat_i[mat_p[p2] + mat_p_count[p2]] = p1 - 1L
                        }
                    } else {
                        mat_p_count[p1] = mat_p_count[p1] + 1L

                        fid_index = mat_p[p1] + mat_p_count[p1]
                        mat_x[fid_index] = 1 - fid[c]
                        mat_i[fid_index] = p1 - 1L
                    }
                }
                row_indices = row_indices[!is.na(row_indices)]
                mat_x[row_indices] = mat_x[row_indices] / rs
            }
        }
    }

    .csc(mat_p, mat_i, mat_x, ncells)
}


#' CRW function
#'
#' An internal function for creating RW samc objects
#'
#' @noRd
.crw = function(x, absorption, fidelity, model) {
    dir = model$dir
    fun = model$fun

    tr = .tr_vals(x, fun, dir)

    kappa_mode = model$kappa_mode

    edge_counts = sum(is.finite(tr))
    # Integer edge numbers, NA for non-edges
    edge_nums = rep(NA_integer_, length(tr))
    edge_nums[is.finite(tr)] = seq_len(edge_counts)


    nrows = terra::nrow(x)
    ncols = terra::ncol(x)
    cell_nums = terra::cells(x)

    # Angle matrix. Rows/cols are in 8-direction order: NW, N, NE, W, E, SW, S, SE
    heading_vecs = matrix(c(-1, 1,
                             0, 1,
                             1, 1,
                            -1, 0,
                             1, 0,
                            -1, -1,
                             0, -1,
                             1, -1),
        nrow = 8, byrow = TRUE)

    # ang_mat was original general version. cos_mat is for von mises, which has a cos() that cancels the acos()
    # ang_mat = matrix(nrow = 8, ncol = 8)
    cos_mat = matrix(nrow = 8, ncol = 8)

    for (r in 1:8) {
        for (c in 1:8) {
            mag_v1 = sqrt(sum(heading_vecs[r, ]^2))
            mag_v2 = sqrt(sum(heading_vecs[c, ]^2))

            # ang_mat[r, c] = acos(sum(heading_vecs[r, ] * heading_vecs[c, ]) / (mag_v1 * mag_v2))
            cos_mat[r, c] = sum(heading_vecs[r, ] * heading_vecs[c, ]) / (mag_v1 * mag_v2)
        }
    }

    # With dir = 4, heading indices 1-4 are N, W, E, S, so keep only those rows/cols
    if (dir == 4) {
        cos_mat = cos_mat[c(2, 4, 5, 7), c(2, 4, 5, 7)]
    }


    # tr offset vars
    dir_vec = .dir_order(dir)
    nb_dirs = dir_vec[dir_vec != 0]
    offsets = .nb_offsets(x, dir) * dir


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
                        p2 = p1 + offsets[c, d]

                        # loop through valid secondary edges
                        for (dv in dir_vec) {
                            if (dv) {
                                e2 = p2 + dv
                                if (!is.na(tr[e2])) {
                                    mat_p[edge_nums[e2] + 1L] = mat_p[edge_nums[e2] + 1L] + 1L
                                }
                            } else {
                                mat_p[edge_nums[e1] + 1L] = mat_p[edge_nums[e1] + 1L] + 1L
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
    crw_index = 0L
    index = 0L

    kappa_vals = terra::values(x, mat = FALSE, row = 1, nrows = 1) * 0

    if (kappa_mode == "scalar") {
        kappa_vals = kappa_vals + model$kappa
    }

    # Loop through cells with values
    for (r in 1:nrows) {
        fid = terra::values(fidelity, mat = FALSE, row = r, nrows = 1)
        vals = 1 - fid - terra::values(absorption, mat = FALSE, row = r, nrows = 1)

        if (kappa_mode == "raster"){
            kappa_vals = terra::values(model$kappa, mat = FALSE, row = r, nrows = 1)
        }

        for (c in 1:ncols) {
            cell = cell + 1
            if (is.finite(vals[c])) {
                p1 = (cell - 1) * dir
                index = index + 1L

                # loop through valid edges
                for (d in 1:dir) {
                    e1 = p1 + d
                    if (!is.na(tr[e1])) {
                        p2 = p1 + offsets[c, d]
                        e1_num = edge_nums[e1]

                        crw_index = crw_index + 1L
                        crw_map[crw_index, ] = c(index, d)

                        rs = 0

                        row_indices = numeric(dir)
                        row_indices[] = NA

                        # Von Mises weights are exp(kappa * cos(turn)) up to a constant that cancels when the row is
                        # normalized by rs. Shifting by the smallest available turn makes the largest weight exactly
                        # 1, so rows can't underflow to 0/0 when straight ahead is unavailable (e.g., at the edge of
                        # the map), and very large kappa converges to the smallest available turn.
                        valid_dv = nb_dirs[!is.na(tr[p2 + nb_dirs])]
                        cos_max = if (length(valid_dv)) max(cos_mat[d, valid_dv]) else 1

                        # loop through valid secondary edges
                        for (dv in dir_vec) {
                            if (dv) {
                                e2 = p2 + dv

                                if (!is.na(tr[e2])) {
                                    e2_num = edge_nums[e2]
                                    mat_p_count[e2_num] = mat_p_count[e2_num] + 1L

                                    temp = exp(kappa_vals[c] * (cos_mat[d, dv] - cos_max))

                                    res = tr[e2] * temp

                                    rs = rs + res

                                    row_indices[dv] = mat_p[e2_num] + mat_p_count[e2_num]
                                    mat_x[mat_p[e2_num] + mat_p_count[e2_num]] = -res * vals[c]
                                    mat_i[mat_p[e2_num] + mat_p_count[e2_num]] = e1_num - 1L
                                }
                            } else {
                                mat_p_count[e1_num] = mat_p_count[e1_num] + 1L

                                fid_index = mat_p[e1_num] + mat_p_count[e1_num]
                                mat_x[fid_index] = 1 - fid[c]
                                mat_i[fid_index] = e1_num - 1L
                            }
                        }
                        row_indices = row_indices[!is.na(row_indices)]
                        mat_x[row_indices] = mat_x[row_indices] / rs
                    }
                }
            }
        }
    }

    mat = .csc(mat_p, mat_i, mat_x, edge_counts)

    dim(tr) = c(dir, length(tr) / dir)
    tr = scale(tr, FALSE, colSums(tr, na.rm = TRUE))
    attr(tr, "scaled:scale") = NULL

    return(
        list(tr = mat,
            crw = crw_map,
            prob = tr,
            abs = terra::values(absorption)[cell_nums[crw_map[, 1]]])
    )
}


#' Transition values
#'
#' Calculate the unnormalized transition value for every cell and direction.
#' The result is a vector with `dir` values per raster cell, in cell order.
#'
#' For RW/CRW models (`ssc = NULL`) the function value is divided by the
#' distance between cells. For SSF models it is multiplied by the distance
#' raised to the power `ssc`.
#'
#' @param data A SpatRaster
#' @param fun A function, or the name of a built-in function: `"1/mean(x)"`
#' (RW/CRW) or `"x[2]"` (SSF)
#' @param dir Either 4 or 8
#' @param ssc Step selection coefficient for SSF models, otherwise `NULL`
#' @noRd
.tr_vals = function(data, fun, dir, ssc = NULL) {
    ssf = !is.null(ssc)

    nrows = terra::nrow(data)
    ncols = terra::ncol(data)

    result = numeric(nrows * ncols * dir)
    index = 0

    dist = .build_lookup_mat(data, dir)
    if (ssf) dist = dist^(ssc)

    if (dir == 4) {
        dir = c(2, 4, 6, 8)
    } else if (dir == 8) {
        dir = c(1:4, 6:9)
    }

    if (is(fun, "function")) {
        for (r in 1:nrows) {
            vals = terra::focalValues(data, 3, r, 1)

            for (c in 1:ncols) {
                v = vals[c, 5]
                for (d in dir) {
                    index = index + 1

                    if (ssf) {
                        result[index] = fun(c(v, vals[c, d])) * dist[r, d]
                    } else {
                        result[index] = fun(c(v, vals[c, d])) / dist[r, d]
                    }
                }
            }
        }
    } else if (!ssf && fun == "1/mean(x)") {
        for (r in 1:nrows) {
            vals = terra::focalValues(data, 3, r, 1)

            for (c in 1:ncols) {
                v = vals[c, 5]
                for (d in dir) {
                    index = index + 1

                    result[index] = 2 / ((v + vals[c, d]) * dist[r, d])
                }
            }
        }
    } else if (ssf && fun == "x[2]") {
        for (r in 1:nrows) {
            vals = terra::focalValues(data, 3, r, 1)

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


.build_lookup_mat = function(rast, dir) {

    # TODO Create tests to directly validate results for both directions options for planar and latlon

    lonlat = terra::is.lonlat(rast)
    nrows = terra::nrow(rast)
    ncols = terra::ncol(rast)

    mat = matrix(nrow = nrows, ncol = 9)
    mat[, c(2, 4, 6, 8)] = 1

    if (dir == 8) {
        mat[, c(1, 3, 7, 9)] = sqrt(2)
    }


    if (!is.na(lonlat)) {
        if (lonlat) {
            cn = (0:(nrows - 1)) * ncols + 1

            adj = terra::adjacent(rast, cn, directions = dir, pairs = TRUE)

            dist = terra::distance(
                terra::xyFromCell(rast, adj[, 1]),
                terra::xyFromCell(rast, adj[, 2]),
                lonlat = TRUE, pairwise = TRUE)

            adj = terra::adjacent(rast, cn, directions = dir, pairs = FALSE)
            adj = t(adj)

            adj[!is.nan(adj)] = dist
            adj = t(adj)

            if (dir == 4) {
                mat[, c(2, 4, 6, 8)] = adj
                mat[, 4] = mat[, 6]
            } else {
                mat[, c(1:4, 6:9)] = adj
                mat[, 1] = mat[, 3]
                mat[, 4] = mat[, 6]
                mat[, 7] = mat[, 9]
            }
        }
    }

    return(mat)
}


#' Convolution
#'
#' TODO description here
#'
#' @param x A function
#' @noRd

.convolution = function(res, abso, fid, dir, sym, threads, precision) {

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

    if (precision == "single") {
        return(.build_convolution_cache_float(kernel, res, fid, abso, sym, threads))
    } else if (precision == "double") {
        return(.build_convolution_cache_double(kernel, res, fid, abso, sym, threads))
    } else {
        stop("Invalid data type. Must be either 'single' or 'double'", call. = FALSE)
    }
}
