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
    nb = .nb_steps(dir)

    offsets = matrix(nb$dr * ncols + nb$dc, nrow = ncols, ncol = dir, byrow = TRUE)

    if (isTRUE(terra::is.lonlat(x, global = TRUE))) {
        offsets[1, nb$dc == -1] = offsets[1, nb$dc == -1] + ncols
        offsets[ncols, nb$dc == 1] = offsets[ncols, nb$dc == 1] - ncols
    }

    offsets
}


#' Neighbor steps
#'
#' Row and column steps (`dr`, `dc`, integer vectors of length `dir`) from a
#' cell to each of its neighbors, in the direction order used by `.tr_vals()`:
#' NW, N, NE, W, E, SW, S, SE (or N, W, E, S for `dir = 4`). This order is
#' ascending in cell order.
#'
#' @noRd
.nb_steps = function(dir) {
    if (dir == 4) {
        list(dr = c(-1L, 0L, 0L, 1L), dc = c(0L, -1L, 1L, 0L))
    } else if (dir == 8) {
        list(dr = c(-1L, -1L, -1L, 0L, 0L, 1L, 1L, 1L), dc = c(-1L, 0L, 1L, -1L, 1L, -1L, 0L, 1L))
    } else {
        stop("Invalid directions", call. = FALSE)
    }
}


#' Strips of raster rows
#'
#' Split `nrows` raster rows into strips of `block_rows` rows (the last one may
#' be shorter; `block_rows` is clamped to `1:nrows`). Returns an integer matrix
#' with one row per strip: the first and last row of the strip (`r1`, `r2`) and
#' of the strip plus a one row halo above and below where those rows exist
#' (`ra`, `rb`).
#'
#' @noRd
.strips = function(nrows, block_rows) {
    nrows = as.integer(nrows)
    block_rows = as.integer(max(1, min(block_rows, nrows)))
    r1 = seq.int(1L, nrows, by = block_rows)
    r2 = pmin(r1 + block_rows - 1L, nrows)

    cbind(r1 = r1, r2 = r2, ra = pmax(1L, r1 - 1L), rb = pmin(nrows, r2 + 1L))
}


#' Neighbor cells within a slab
#'
#' For destination cells at slab rows `drow` and columns `dcol`, the slab cell
#' index (row-major, 1-based) of the neighbor one step of `dr` rows and `dc`
#' columns away, or NA where that neighbor is outside the slab. With `wrap`,
#' columns wrap around (global lonlat rasters), so only rows can fall outside.
#'
#' @noRd
.slab_nb = function(drow, dcol, dr, dc, nh, ncols, wrap) {
    qr = drow + dr
    qc = dcol + dc
    if (wrap) qc = (qc - 1L) %% ncols + 1L

    q = (qr - 1L) * ncols + qc
    q[qr < 1L | qr > nh | qc < 1L | qc > ncols] = NA
    q
}


#' Sort entry slots by row index
#'
#' `slot_x` (values) and `slot_i` (row indices) are slot x column matrices for
#' some matrix columns, NA in `slot_i` marking an empty slot. Returns the values
#' and indices with each column's slots sorted by row index, empty slots last,
#' as vectors in the same column-major layout.
#'
#' @noRd
.sort_slots = function(slot_x, slot_i) {
    o = order(rep(seq_len(ncol(slot_i)), each = nrow(slot_i)), slot_i)
    list(x = slot_x[o], i = slot_i[o])
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
#' An internal function for creating the matrix for RW and SSF samc objects.
#' States are the non-NA cells of the map, numbered in cell order.
#'
#' The matrix is built in strips of whole raster rows. This function does all
#' of the raster access and bookkeeping; each strip is handed to `.rw_block()`
#' as plain vectors, and the finished slab of the CSC arrays it returns is
#' copied into the preallocated result.
#'
#' @param x A SpatRaster with the data for the transition function
#' @param absorption A SpatRaster with the total absorption
#' @param fidelity A SpatRaster with the fidelity
#' @param model The model list. Uses `name`, `fun`, `dir`, and for SSF, `ssc`
#' @param block_cells Approximate number of destination cells per strip (whole
#' rows, at least one). The strip temporaries are about 1 kB per cell at
#' `dir = 8`.
#'
#' @return A dgCMatrix
#' @noRd
.rw = function(x, absorption, fidelity, model, block_cells = 2e4) {
    dir = as.integer(model$dir)

    if (model$name == "RW") {
        tr = .tr_vals(x, model$fun, dir)
    } else if (model$name == "SSF") {
        tr = .tr_vals(x, model$fun, dir, ssc = model$ssc)
    } else {
        stop("Unexpected model in .rw()", call. = FALSE)
    }

    nrows = as.integer(terra::nrow(x))
    ncols = as.integer(terra::ncol(x))
    row_tr = as.numeric(ncols) * dir # tr values per raster row (double, so that tr indices can't overflow)
    nb = .nb_steps(dir)
    wrap = isTRUE(terra::is.lonlat(x, global = TRUE))
    strips = .strips(nrows, block_cells %/% ncols)

    # State numbers (1-based in cell order, 0 for NA cells), and the number of non-zeros: a diagonal entry per state
    # plus an entry per transition value out of a state. The count is exact for the built-in transition functions,
    # whose values are NA wherever the destination is NA or off the map. A blocked scan rather than terra::cells(),
    # which reads the whole raster at once.
    state = integer(nrows * ncols)
    nstates = 0L
    nnz = 0
    for (s in seq_len(nrow(strips))) {
        r1 = strips[s, "r1"]
        r2 = strips[s, "r2"]

        valid = !is.na(terra::values(x, mat = FALSE, row = r1, nrows = r2 - r1 + 1L))
        k = sum(valid)
        if (k > 0L) {
            state[((r1 - 1L) * ncols + 1L):(r2 * ncols)][valid] = nstates + seq_len(k)

            # Transition values of the block, less those out of NA cells
            t0 = (r1 - 1) * row_tr
            na_tr = t0 + rep((which(!valid) - 1) * dir, each = dir) + seq_len(dir)
            nnz = nnz + k + sum(!is.na(tr[(t0 + 1):(r2 * row_tr)])) - sum(!is.na(tr[na_tr]))
        }
        nstates = nstates + k
    }
    rm(valid)
    if (nnz > .Machine$integer.max) stop("Too many non-zero entries for a dgCMatrix", call. = FALSE)

    invisible(gc(full = FALSE)) # drop the scan garbage before the large allocations

    mat_x = numeric(nnz)
    mat_i = integer(nnz)
    mat_p = integer(nstates + 1L)
    pos = 0L # entries written so far
    col_pos = 0L # columns (states) written so far

    for (s in seq_len(nrow(strips))) {
        r1 = strips[s, "r1"]
        r2 = strips[s, "r2"]
        ra = strips[s, "ra"]
        rb = strips[s, "rb"]

        blk = .rw_block(
            tr = tr[((ra - 1) * row_tr + 1):(rb * row_tr)],
            fid = terra::values(fidelity, mat = FALSE, row = ra, nrows = rb - ra + 1L),
            abs = terra::values(absorption, mat = FALSE, row = ra, nrows = rb - ra + 1L),
            state = state[((ra - 1L) * ncols + 1L):(rb * ncols)],
            ncols = ncols,
            dest_rows = c(r1, r2) - ra + 1L,
            dr = nb$dr,
            dc = nb$dc,
            wrap = wrap)

        m = length(blk$x)
        if (m > 0L) {
            idx = pos + seq_len(m)
            mat_x[idx] = blk$x
            mat_i[idx] = blk$i
        }
        ns = length(blk$p)
        mat_p[col_pos + 1L + seq_len(ns)] = pos + cumsum(blk$p)
        pos = pos + m
        col_pos = col_pos + ns
        rm(blk)

        invisible(gc(full = FALSE))
    }

    # Only a custom transition function with values towards NA cells can leave unused slots (see the count above)
    if (pos < nnz) {
        length(mat_x) = pos
        length(mat_i) = pos
    }

    .csc(mat_p, mat_i, mat_x, nstates)
}


#' RW matrix block
#'
#' Build one slab of the RW/SSF matrix: the CSC columns of the destination
#' cells of a strip. Takes only plain data and does no raster access, so that
#' it can be ported to C++ or used on standalone tiles.
#'
#' The caller supplies a slab of `nh` rows of `ncols` cells, row-major (cell
#' order). The destinations are the cells of slab rows
#' `dest_rows[1]:dest_rows[2]`. Entries come from the neighbors of the
#' destinations, so the slab must include a one row halo: every map row within
#' one row of the destination rows. Slab rows beyond the first and last are
#' treated as off the map, and so are columns beyond the first and last unless
#' `wrap` is set.
#'
#' @param tr Transition values of the slab cells, `dir` per cell, in cell order
#' and the direction order of `dr`/`dc` (as from `.tr_vals()`). NA where there
#' is no transition.
#' @param fid,abs Fidelity and absorption of the slab cells
#' @param state Integer state number of each slab cell, 0 if the cell is not a
#' state. These are the (1-based) matrix row and column numbers, so they can be
#' global. They must be ascending in cell order (wrapped columns excepted) for
#' the row indices of each column to come out sorted.
#' @param ncols Number of slab columns
#' @param dest_rows First and last destination row, as slab rows
#' @param dr,dc Integer row and column steps from a cell to its neighbors, in
#' ascending cell order and with every step's opposite included (see
#' `.nb_steps()`)
#' @param wrap Whether the first and last slab columns are neighbors (global
#' lonlat raster)
#'
#' @return A list with the slab of the CSC arrays for the destination cells
#' that are states (`state != 0`), in cell order: `p`, the number of entries in
#' each column (integer); `i`, the row indices, 0-based (`state - 1`); and `x`,
#' the values. Within each column the row indices are ascending.
#' @noRd
.rw_block = function(tr, fid, abs, state, ncols, dest_rows, dr, dc, wrap) {
    dir = length(dr)
    nh = length(state) %/% ncols

    # Cells that are not states give no entries
    vals = 1 - fid - abs
    vals[!is.finite(vals) | state == 0L] = NA

    # Row sums of the transition values of each slab cell, added in direction order (colSums() rounds differently)
    rs = numeric(length(state))
    for (j in seq_len(dir)) {
        t = tr[seq.int(j, by = dir, length.out = length(state))]
        t[is.na(t)] = 0
        rs = rs + t
    }

    dest = ((dest_rows[1] - 1L) * ncols + 1L):(dest_rows[2] * ncols)
    n = length(dest)
    drow = rep(dest_rows[1]:dest_rows[2], each = ncols)
    dcol = rep(seq_len(ncols), length.out = n)

    # Gather by destination: the column of a destination has an entry for each neighbor k that is a source, using
    # the neighbor's transition value in the opposite direction, opp[k]. The first h neighbors come before the
    # destination in cell order, so with the diagonal entry in slot h + 1 the dir + 1 slots of a column are in
    # ascending row order (except for wrapped columns, sorted below).
    opp = match(paste(-dr, -dc), paste(dr, dc))
    h = sum(dr < 0L | (dr == 0L & dc < 0L))

    slot_x = matrix(NA_real_, dir + 1L, n) # values, one column per destination
    slot_i = matrix(NA_integer_, dir + 1L, n) # row numbers (source states), NA for no entry
    for (k in seq_len(dir)) {
        q = .slab_nb(drow, dcol, dr[k], dc[k], nh, ncols, wrap)
        tq = tr[(q - 1L) * dir + opp[k]]
        iq = state[q]
        iq[is.na(tq) | is.na(vals[q])] = NA

        slot = k + (k > h)
        slot_x[slot, ] = (-tq * vals[q]) / rs[q]
        slot_i[slot, ] = iq
    }
    slot_x[h + 1L, ] = 1 - fid[dest]
    slot_i[h + 1L, ] = state[dest]

    # On a global raster the W/E neighbors of the first/last column are on the other side of the map, out of cell
    # order
    if (wrap) {
        wc = which(dcol == 1L | dcol == ncols)
        s = .sort_slots(slot_x[, wc, drop = FALSE], slot_i[, wc, drop = FALSE])
        slot_x[, wc] = s$x
        slot_i[, wc] = s$i
        rm(s)
    }

    is_col = state[dest] != 0L
    keep = !is.na(slot_i) & rep(is_col, each = dir + 1L)
    p = as.integer(colSums(keep))[is_col]
    x = slot_x[keep]
    rm(slot_x)
    i = slot_i[keep] - 1L

    list(p = p, i = i, x = x)
}


#' CRW function
#'
#' An internal function for creating the matrix for CRW samc objects. States
#' are the edges of the map (a cell and the heading it is left with) rather
#' than the cells: the (cell, direction) pairs with a finite transition value,
#' numbered in that order.
#'
#' The matrix is built in strips of whole raster rows, as in `.rw()`, with each
#' strip handed to `.crw_block()`.
#'
#' @param x A SpatRaster with the data for the transition function
#' @param absorption A SpatRaster with the total absorption
#' @param fidelity A SpatRaster with the fidelity
#' @param model The model list. Uses `fun`, `dir`, `kappa`, and `kappa_mode`
#' @param block_entries Approximate number of entry slots per strip (whole
#' rows, at least one). A strip has `(dir + 1) * dir` slots per cell.
#'
#' @return A list with the matrix (`tr`, a dgCMatrix), a two column matrix
#' with the cell and heading of each state (`crw`), the normalized transition
#' values as a `dir` x cells matrix (`prob`), and the absorption for each state
#' (`abs`)
#' @noRd
.crw = function(x, absorption, fidelity, model, block_entries = 2.5e5) {
    dir = as.integer(model$dir)
    tr = .tr_vals(x, model$fun, dir)

    nrows = as.integer(terra::nrow(x))
    ncols = as.integer(terra::ncol(x))
    ncell = nrows * ncols
    row_tr = as.numeric(ncols) * dir # tr values per raster row (double, so that tr indices can't overflow)
    nb = .nb_steps(dir)
    wrap = isTRUE(terra::is.lonlat(x, global = TRUE))
    kappa_raster = model$kappa_mode == "raster"

    # Cosine of the turn angle between each pair of headings, in direction order. The heading of a step is (dc, -dr),
    # i.e. x to the right and y up. The von Mises weights only need the cosine, so the angle itself (the acos() of
    # this) is never calculated.
    heading = cbind(nb$dc, -nb$dr)
    mag = sqrt(rowSums(heading^2))
    cos_mat = tcrossprod(heading) / tcrossprod(mag)

    # Cells with entries (finite 1 - fidelity - absorption), and the cell and direction of each edge
    abs_all = terra::values(absorption, mat = FALSE)
    has_entries = is.finite(1 - terra::values(fidelity, mat = FALSE) - abs_all)
    e = which(is.finite(tr))
    nedges = length(e)
    cell = (e - 1L) %/% dir + 1L
    d = e - (cell - 1L) * dir
    rm(e)

    # Edges before each raster row, so that a strip can number its edges from its own slab of tr
    row_base = c(0L, cumsum(tabulate((cell - 1L) %/% ncols + 1L, nrows)))

    # Number of non-zeros: column (p2, dv) has an entry for each edge into p2 from a cell with entries, plus the
    # diagonal if p2 has entries itself. The whole map is treated as one slab to find where each edge goes.
    from = has_entries[cell]
    src = cell[from]
    dd = d[from]
    dest = .slab_nb((src - 1L) %/% ncols + 1L, (src - 1L) %% ncols + 1L, nb$dr[dd], nb$dc[dd], nrows, ncols, wrap)
    indeg = tabulate(dest, ncell)
    rm(src, dd, dest)
    nnz = sum(as.numeric(tabulate(cell, ncell)) * (indeg + has_entries))
    rm(indeg)
    if (nnz > .Machine$integer.max) stop("Too many non-zero entries for a dgCMatrix", call. = FALSE)

    # The cell (as its index among cells with entries) and direction of each state: filled for the edges of cells
    # with entries, in order, and left as 0 in the remaining rows at the end
    crw_map = matrix(0L, nrow = nedges, ncol = 2)
    m = sum(from)
    crw_map[seq_len(m), 1] = cumsum(has_entries)[cell[from]]
    crw_map[seq_len(m), 2] = d[from]
    rm(cell, d, from, has_entries)
    abs_out = abs_all[terra::cells(x)[crw_map[, 1]]]
    rm(abs_all)

    mat_x = numeric(nnz)
    mat_i = integer(nnz)
    mat_p = integer(nedges + 1L)
    pos = 0L # entries written so far
    col_pos = 0L # columns (states) written so far

    strips = .strips(nrows, block_entries %/% ((dir + 1L) * dir * ncols))
    for (s in seq_len(nrow(strips))) {
        r1 = strips[s, "r1"]
        r2 = strips[s, "r2"]
        ra = strips[s, "ra"]
        rb = strips[s, "rb"]

        # Edge numbers of the slab cells: NA for non-edges, otherwise counting on from the edges of earlier rows
        ts = tr[((ra - 1) * row_tr + 1):(rb * row_tr)]
        ev = is.finite(ts)
        edge = rep(NA_integer_, length(ts))
        edge[ev] = row_base[ra] + seq_len(sum(ev))
        rm(ev)

        kappa = model$kappa
        if (kappa_raster) kappa = terra::values(kappa, mat = FALSE, row = ra, nrows = rb - ra + 1L)

        blk = .crw_block(
            tr = ts,
            edge = edge,
            fid = terra::values(fidelity, mat = FALSE, row = ra, nrows = rb - ra + 1L),
            abs = terra::values(absorption, mat = FALSE, row = ra, nrows = rb - ra + 1L),
            kappa = kappa,
            cos_mat = cos_mat,
            ncols = ncols,
            dest_rows = c(r1, r2) - ra + 1L,
            dr = nb$dr,
            dc = nb$dc,
            wrap = wrap)
        rm(ts, edge)

        m = length(blk$x)
        if (m > 0L) {
            idx = pos + seq_len(m)
            mat_x[idx] = blk$x
            mat_i[idx] = blk$i
        }
        ns = length(blk$p)
        mat_p[col_pos + 1L + seq_len(ns)] = pos + cumsum(blk$p)
        pos = pos + m
        col_pos = col_pos + ns
        rm(blk)

        invisible(gc(full = FALSE)) # see .rw()
    }

    if (pos != nnz || col_pos != nedges) {
        stop("Unexpected number of matrix entries in .crw(). Please report with a minimum reproducible example.",
            call. = FALSE)
    }

    # Normalized transition values. Computed after the strips, which gave a lower memory peak.
    dim(tr) = c(dir, ncell)
    tr = tr / rep(colSums(tr, na.rm = TRUE), each = dir)

    list(tr = .csc(mat_p, mat_i, mat_x, nedges), crw = crw_map, prob = tr, abs = abs_out)
}


#' CRW matrix block
#'
#' Build one slab of the CRW matrix: the CSC columns of the edges (states)
#' leaving the destination cells of a strip. Takes only plain data and does no
#' raster access, so that it can be ported to C++ or used on standalone tiles.
#'
#' The slab and halo requirements are the same as for `.rw_block()`.
#'
#' The entry from row edge (q, dd), which enters destination cell p2 from
#' neighbor q, to column edge (p2, dv) is
#' `-tr[p2, dv] * w(dd, dv) * (1 - fid[q] - abs[q]) / rs`, where `w` is the
#' von Mises turn weight with the kappa of q and `rs` normalizes over the edges
#' leaving p2. The diagonal entry of column (p2, dv) is `1 - fid[p2]`.
#'
#' @param tr Transition values of the slab cells, `dir` per cell, in cell order
#' and the direction order of `dr`/`dc` (as from `.tr_vals()`)
#' @param edge Integer edge (state) number of each slab cell and direction,
#' laid out as `tr`, NA where that is not an edge. These are the (1-based)
#' matrix row and column numbers, so they can be global. They must be ascending
#' in cell order (wrapped columns excepted) for the row indices of each column
#' to come out sorted.
#' @param fid,abs Fidelity and absorption of the slab cells
#' @param kappa A single kappa, or the kappa of each slab cell
#' @param cos_mat `dir` x `dir` matrix of the cosine of the turn angle between
#' each pair of directions
#' @param ncols,dest_rows,dr,dc,wrap As for `.rw_block()`
#'
#' @return A list with the slab of the CSC arrays for the edges of the
#' destination cells, in edge order: `p`, the number of entries in each column
#' (integer); `i`, the row indices, 0-based (`edge - 1`); and `x`, the values.
#' Within each column the row indices are ascending.
#' @noRd
.crw_block = function(tr, edge, fid, abs, kappa, cos_mat, ncols, dest_rows, dr, dc, wrap) {
    dir = length(dr)
    nh = length(fid) %/% ncols

    # Cells without a finite value give no entries
    vals = 1 - fid - abs
    vals[!is.finite(vals)] = NA

    dest = ((dest_rows[1] - 1L) * ncols + 1L):(dest_rows[2] * ncols)
    n = length(dest)
    drow = rep(dest_rows[1]:dest_rows[2], each = ncols)
    dcol = rep(seq_len(ncols), length.out = n)

    # Neighbor k of a destination reaches it in direction opp[k]; the first h neighbors come before it in cell
    # order (see .rw_block())
    opp = match(paste(-dr, -dc), paste(dr, dc))
    h = sum(dr < 0L | (dr == 0L & dc < 0L))

    # The edges leaving the destination cells are the matrix columns. e2 is their numbers in column order (direction
    # within destination cell), NA where (p2, dv) is not an edge; t2 is their transition values as a destination x
    # direction matrix, NA where there is no edge.
    di = ((dest[1] - 1L) * dir + 1L):(dest[n] * dir)
    e2 = edge[di]
    t2 = tr[di]
    t2[is.na(e2)] = NA
    dim(t2) = c(dir, n)
    t2 = t(t2)

    # One matrix column per (destination, dv), one row per entry slot as in .rw_block(): with no column wrap the
    # slots are ascending in row (edge) number, since edges are numbered in cell order
    slot_x = matrix(NA_real_, dir + 1L, dir * n)
    slot_i = matrix(NA_integer_, dir + 1L, dir * n)
    for (k in seq_len(dir)) {
        dd = opp[k] # heading of the row edge (q, dd) from neighbor q into the destination
        q = .slab_nb(drow, dcol, dr[k], dc[k], nh, ncols, wrap)
        eq = edge[(q - 1L) * dir + dd]
        eq[is.na(vals[q])] = NA
        kq = if (length(kappa) == 1L) kappa else kappa[q]

        # Von Mises weights are exp(kappa * cos(turn)) up to a constant that cancels when the row is normalized by
        # rs. Shifting by the smallest available turn (the largest cosine over the edges leaving p2, 1 if there are
        # none) makes the largest weight exactly 1, so rows can't underflow to 0/0 when straight ahead is unavailable
        # (e.g., at the edge of the map), and very large kappa converges to the smallest available turn.
        cm = cos_mat[dd, ]
        cmax = rep(-Inf, n)
        for (dv in seq_len(dir)) {
            u = !is.na(t2[, dv]) & cm[dv] > cmax
            cmax[u] = cm[dv]
        }
        cmax[cmax == -Inf] = 1

        # tr * exp(kappa * (cos - cos_max)), the row sum accumulated in dv order, then (-w * vals) / rs. Changing
        # the order of operations would change the values in the last bits.
        wt = t2 * exp(kq * (rep(cm, each = n) - cmax))
        rs = numeric(n)
        for (dv in seq_len(dir)) {
            w = wt[, dv]
            w[is.na(w)] = 0
            rs = rs + w
        }

        ie = rep(eq, each = dir)
        ie[is.na(e2)] = NA

        slot = k + (k > h)
        slot_x[slot, ] = t((-wt * vals[q]) / rs)
        slot_i[slot, ] = ie
    }
    rm(wt, w, t2)

    ie = e2
    ie[rep(is.na(vals[dest]), each = dir)] = NA
    slot_x[h + 1L, ] = rep(1 - fid[dest], each = dir)
    slot_i[h + 1L, ] = ie

    # On a global raster the W/E neighbors of the first/last column are on the other side of the map, out of edge
    # order
    if (wrap) {
        wcell = which(dcol == 1L | dcol == ncols)
        wc = rep((wcell - 1L) * dir, each = dir) + seq_len(dir)
        s = .sort_slots(slot_x[, wc, drop = FALSE], slot_i[, wc, drop = FALSE])
        slot_x[, wc] = s$x
        slot_i[, wc] = s$i
        rm(s)
    }

    keep = !is.na(slot_i)
    p = as.integer(colSums(keep))[!is.na(e2)]
    x = slot_x[keep]
    rm(slot_x)
    i = slot_i[keep] - 1L

    list(p = p, i = i, x = x)
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
#' @param block_cells Approximate number of cells read per
#' `terra::focalValues()` call. Each block is processed one direction at a
#' time, so this bounds the temporary memory (roughly 110 bytes per cell of the
#' block) without paying the per-call overhead of reading one row at a time.
#' @noRd
.tr_vals = function(data, fun, dir, ssc = NULL, block_cells = 2.5e5) {
    ssf = !is.null(ssc)

    if (is.function(fun)) {
        mode = "custom"
    } else if (!ssf && fun == "1/mean(x)") {
        mode = "rw"
    } else if (ssf && fun == "x[2]") {
        mode = "ssf"
    } else {
        stop("Invalid transition function defined", call. = FALSE)
    }

    nrows = terra::nrow(data)
    ncols = terra::ncol(data)

    dist = .build_lookup_mat(data, dir)
    if (ssf) dist = dist^(ssc)

    # focalValues() columns of the neighbors, in direction order
    nb_cols = if (dir == 4) c(2, 4, 6, 8) else c(1:4, 6:9)

    # One column per cell and one row per direction, so that flattening it at the end puts the directions of each
    # cell next to each other
    result = matrix(0, nrow = dir, ncol = nrows * ncols)

    block_rows = max(1L, as.integer(block_cells %/% ncols))

    for (r1 in seq(1L, nrows, by = block_rows)) {
        rows = r1:min(r1 + block_rows - 1L, nrows)
        cells = ((r1 - 1L) * ncols + 1L):(rows[length(rows)] * ncols)
        cell_rows = rep(rows, each = ncols)

        # One row per cell of the block. Column 5 is the cell itself, nb_cols are its neighbors (NA if off the map)
        vals = terra::focalValues(data, 3, r1, length(rows))
        from = vals[, 5]

        for (j in seq_len(dir)) {
            to = vals[, nb_cols[j]]
            step = dist[cell_rows, nb_cols[j]]

            result[j, cells] = switch(mode,
                rw = 2 / ((from + to) * step),
                ssf = to * step,
                custom = {
                    # A plain loop rather than vapply(), which would add a second function call per pair
                    v = numeric(length(from))
                    for (k in seq_along(from)) v[k] = fun(c(from[k], to[k]))
                    if (ssf) v * step else v / step
                })
        }
    }

    dim(result) = NULL
    result
}


#' Build distance lookup matrix
#'
#' Distances from a cell to each of its neighbors. Columns are the positions
#' of a 3x3 window in the same order as `terra::focalValues()`. The center
#' column and, for `dir = 4`, the diagonal columns are NA. For lonlat rasters
#' the distances depend on latitude, so there is one row per raster row.
#'
#' @param rast A SpatRaster
#' @param dir Either 4 or 8
#'
#' @return A matrix with `nrow(rast)` rows and 9 columns
#' @noRd
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


    if (isTRUE(lonlat)) {
        # Distances only vary by row, so use the first cell of each row
        cn = (0:(nrows - 1)) * ncols + 1

        adj = terra::adjacent(rast, cn, directions = dir, pairs = TRUE)

        dist = terra::distance(
            terra::xyFromCell(rast, adj[, 1]),
            terra::xyFromCell(rast, adj[, 2]),
            lonlat = TRUE, pairwise = TRUE)

        # `dist` has one value per existing neighbor, ordered by cell and then by direction. With `pairs = FALSE` the
        # result is a cell x direction matrix with NaN for missing neighbors. R fills matrices by column, so
        # transpose to direction x cell to fill the non-NaN positions in the same order as `dist`, then transpose back.
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

    return(mat)
}


#' Convolution
#'
#' An internal function for creating the cache used by the convolution
#' algorithm. Non-finite values in the inputs are treated as 0.
#'
#' @param res A SpatRaster with the resistance data
#' @param abso A SpatRaster with the total absorption
#' @param fid A SpatRaster with the fidelity
#' @param dir Either 4 or 8
#' @param sym Whether the transition function is symmetric
#' @param threads Number of threads to use
#' @param precision Either `"single"` or `"double"`
#'
#' @return The convolution cache created by the C++ code
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

    nr = terra::nrow(res)
    nc = terra::ncol(res)

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
