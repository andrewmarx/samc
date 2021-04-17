# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R check.R
NULL


#' Create an samc object
#'
#' Create an samc object that contains the absorbing Markov chain data
#'
#' This function is used to create a \code{\link{samc-class}} object. There are
#' two options for creating this object.
#'
#' \strong{Option 1: Raster and Matrix Inputs}
#'
#' The \code{\link{samc-class}} object can be created from a combination of
#' resistance, absorption, and fidelity data. These different landscape data
#' inputs must be the same type (a matrix or RasterLayer), and have identical
#' properties, including dimensions, location of NA cells, and CRS (if using
#' RasterLayers). Some of the inputs are mandatory, whereas others are optional.
#'
#' The \code{resistance} and \code{absorption} inputs are always mandatory, whereas the
#' \code{fidelity} input is optional. If the \code{fidelity} input is not provided, then it it
#' is assumed that there is no site fidelity (i.e., individuals will always move
#' to an adjacent cell each time step).
#'
#' The \code{tr_fun} parameter is mandatory. It used when calculating the values for
#' the transition matrix. Internally, this is passed to the \code{\link[gdistance]{transition}}
#' function in the gdistance package to create the transition matrix.
#'
#' \strong{Option 2: P Matrix Input}
#'
#' The \code{p_mat} parameter can be used to create a \code{\link{samc-class}} object
#' directly from a preconstructed P matrix. This matrix must be either a base R
#' matrix, or a sparse matrix (dgCMatrix format) from the Matrix package. It
#' must meet the requirement of a P matrix described in Fletcher et al. (2019).
#' This includes:
#' \itemize{
#'   \item The number of rows must equal the number of columns (a square matrix)
#'   \item The last row must contain all 0's, except the last element, which must be 1
#'   \item Each row must sum to 1
#'   \item All values must be in the range of 0-1
#' }
#'
#' Additionally, the columns and rows of the P matrix may be named (e.g., using
#' dimnames(), rowname(), colnames(), etc). When specifying \code{origin} or \code{dest} inputs
#' to metrics, these names may be used instead of cell numbers. This has the
#' advantage of making the code for an analysis easier to read and interpret,
#' which may also help to eliminate unintentional mistakes. There are two
#' requirements for naming the rows/cols of a P matrix. First, since the P matrix
#' represents a pairwise matrix, the row and column names must be the same. Second,
#' there must be no duplicate names. The exception to these rules is the very last
#' column and the very last row of the P matrix. Since these are not part of the
#' pairwise transition matrix, they may have whatever names the user prefers.
#'
#' \strong{Other Parameters}
#'
#' The \code{directions} parameter is optional. When constructing the P matrix from
#' matrix or raster data, the \code{samc()} function must decide how adjacent cells are
#' connected. This value can be set to either 4 or 8. When set to 4, nodes are
#' connected horizontally and vertically (similar to the directions of how a rook
#' moves in chess). When set to 8, nodes are connected diagonally in addition to
#' horizontally and vertically (queen movement in chess). When not specified,
#' the \code{samc()} function defaults to a value of 8. When using large datasets to
#' construct a P matrix, setting the directions to 4 may lead to significant
#' improvements in computation time and the amount of memory needed to perform
#' an analysis.
#'
#' \strong{Additional Information}
#'
#' Depending on the data used to construct the samc-class object, some metrics
#' may cause crashes. This is a result of the underlying P matrix having specific
#' properties that make some equations unsolvable. One known case is a P matrix
#' that represents a disconnected graph, which can lead to the \code{cond_passage()}
#' function crashing. In terms of raster/matrix inputs, a disconnected graph
#' occurs when one or more pixels/cells are unreachable from other pixels/cells
#' due to the presence of a full barrier made up of NA values. In a raster, these
#' may be obvious as islands, but can be as inconspicuous as a rogue isolated
#' pixel. There may be other currently unknown situations that lead to unsolvable
#' metrics.
#'
#' Future work is planned towards identifying these issues during creation of the
#' samc-class object and handling them appropriately to prevent inadvertent
#' crashes.
#'
#'
#' @param data A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}} or Matrix package dgCMatrix sparse matrix.
#' @param absorption A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param fidelity A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param tr_fun A function to calculate the transition values in the \code{\link[gdistance]{transition}} function
#' @param directions Optional param. Must be set to either 4 or 8 (default is 8)
#' @param symm Optional param for specifying if the transition matrix should be symmetric. Defaults to \code{TRUE}
#'
#' @param resistance Deprecated. Use the `data` parameter.
#' @param p_mat Deprecated. Use the `data` parameter.
#' @param latlon Deprecated. No longer needed.
#' @param override Deprecated. See \code{\link{samc-class}} for the alternative.
#' @param ... Placeholder
#'
#' @return A spatial absorbing Markov chain object
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "samc",
  function(data, absorption, fidelity, tr_fun, ...) {
    standardGeneric("samc")
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_fun = "function"),
  function(data, absorption, fidelity, tr_fun, directions = 8, symm = TRUE, latlon, override) {

    if (!missing(override))
      warning("The override parameter is deprecated. See the samc-class documentation for more details.", call. = FALSE)

    if (!missing(latlon)) {
      warning("latlon is deprecated and no longer needed; please remove it.", call. = FALSE)
    }

    # Make sure the input data all aligns
    check(data, absorption)
    check(data, fidelity)

    if (any(data[] <= 0, na.rm = TRUE)) {
      stop("The data must not have values <= 0", call. = FALSE)
    }

    if (any(absorption[] < 0, na.rm = TRUE)) {
      stop("The absorption data must not have values <= 0", call. = FALSE)
    }

    if (any(absorption[] > 1, na.rm = TRUE)) {
      stop("The absorption data must not have values > 1", call. = FALSE)
    }

    if (sum(absorption[], na.rm = TRUE) == 0) {
      stop("At least one cell must have an absorption value > 0", call. = FALSE)
    }

    if (any(fidelity[] < 0, na.rm = TRUE)) {
      stop("The fidelity data must not have values < 0", call. = FALSE)
    }

    if (any(fidelity[] > 1, na.rm = TRUE)) {
      stop("The fidelity data must not have values > 1", call. = FALSE)
    }

    if (any((fidelity[] + absorption[]) > 1, na.rm = TRUE)) {
      stop("No cells can have fidelity + absoprtion > 1", call. = FALSE)
    }

    if (!(directions %in% c(4, 8))) {
      stop("directions must be set to either 4 or 8", call. = FALSE)
    }

    # Create map template
    m <- data
    m[] <- is.finite(m[])

    # Check for "clumps"
    cl <- raster::clump(m, directions = directions, gaps = FALSE)
    clumps <- sum(!is.na(unique(cl[])))

    if (clumps > 1) {
      print("Warning: Input contains disconnected regions. This does not work with the cond_passage() metric.")

      temp_abs <- absorption
      temp_abs[temp_abs > 0] <- 1
      temp_abs <- temp_abs * cl

      if (!all(1:clumps %in% unique(temp_abs[]))) stop("All disconnected regions must have at least one non-zero absorption value", call. = FALSE)
    }


    # Create the transition matrix
    tr <- gdistance::transition(data, transitionFunction = tr_fun, directions, symm = symm)
    if(directions == 8 || raster::isLonLat(data)) {
      tr <- gdistance::geoCorrection(tr, type = "c")
    }

    if(is.na(raster::projection(data)) && raster::xres(data) != raster::yres(data)) {
      warning("Raster cells are not square (number of columns/rows is not propotional to the spatial extents). There is no defined projection to account for this, so the geocorrection may lead to distortion if the intent was for the raster cells to represent a uniformly spaced grid.", call. = FALSE)
    }



    tr_mat <- gdistance::transitionMatrix(tr)

    # Normalize the transition Matrix
    abs_vec <- as.vector(absorption)
    fid_vec <- as.vector(fidelity)

    tr_mat <- methods::as(tr_mat, "dgTMatrix") # dgTMatrix is easier to edit directly

    Matrix::diag(tr_mat) <- 0
    tr_mat@x <- (1 - abs_vec[tr_mat@i + 1] - fid_vec[tr_mat@i + 1]) * tr_mat@x / Matrix::rowSums(tr_mat)[tr_mat@i + 1]
    Matrix::diag(tr_mat) <- fid_vec

    # Combine the transition matrix with the absorbing data and convert back to dgCmatrix
    samc_df <- data.frame(i = c(tr_mat@i, (0:(length(abs_vec) - 1))[is.finite(abs_vec)], length(abs_vec)),
                         j = c(tr_mat@j, rep(length(abs_vec), sum(is.finite(abs_vec))), length(abs_vec)),
                         x = c(tr_mat@x, abs_vec[is.finite(abs_vec)], 1))

    p = Matrix::sparseMatrix(i = samc_df$i,
                             j = samc_df$j,
                             x = samc_df$x,
                             index1 = FALSE)

    # Adjust fidelity values for isolated cells
    Matrix::diag(p) <- Matrix::diag(p) - Matrix::rowSums(p) + 1

    # Remove rows/cols for NA cells
    excl <- which(is.na(abs_vec))
    if (length(excl) > 0) p = p[-excl, -excl]

    # Check dimnames
    if (is.null(rownames(p))) rownames(p) <- 1:nrow(p)
    if (is.null(colnames(p))) colnames(p) <- 1:ncol(p)

    if (any(duplicated(rownames(p))))
      stop("Row names must be unique")
    if (any(duplicated(colnames(p))))
      stop("Column names must be unique")

    # Assemble final
    samc_mat <- methods::new("samc", p = p, source = "map", map = m, clumps = clumps, override = FALSE)

    return(samc_mat)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_fun = "function"),
  function(data, absorption, tr_fun, directions = 8, symm = TRUE, latlon, override) {

    fidelity <- data
    fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity, tr_fun, directions = directions, symm = symm, latlon = latlon, override = override))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "matrix",
            tr_fun = "function"),
  function(data, absorption, fidelity, tr_fun, directions = 8, symm = TRUE, override) {

    data <- .rasterize(data)
    absorption <- .rasterize(absorption)
    fidelity <- .rasterize(fidelity)

    #fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity, tr_fun, directions = directions, symm = symm, override = override))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "missing",
            tr_fun = "function"),
  function(data, absorption, tr_fun, directions = 8, symm = TRUE, override) {

    data <- .rasterize(data)
    absorption <- .rasterize(absorption)

    return(samc(data, absorption, tr_fun = tr_fun, directions = directions, symm = symm, override = override))
  })

# TODO: stop-gap for parameter changes. Remove in future version.
#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "matrix",
            fidelity = "missing",
            tr_fun = "function"),
  function(absorption, fidelity, tr_fun, directions = 8, symm = TRUE, override, resistance) {

    if (!missing(resistance)) {
      warning("The resistance parameter is depcrecated. Use the data parameter instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data, absorption, tr_fun = tr_fun, directions = directions, symm = symm, override = override))
  })

# TODO: stop-gap for parameter changes. Remove in future version.
#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "RasterLayer",
            fidelity = "missing",
            tr_fun = "function"),
  function(absorption, fidelity, tr_fun, directions = 8, symm = TRUE, override, resistance, latlon) {

    if (!missing(resistance)) {
      warning("The resistance parameter is depcrecated. Use the data parameter instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data, absorption, tr_fun = tr_fun, directions = directions, symm = symm, latlon = latlon, override = override))
  })

# TODO: stop-gap for parameter changes. Remove in future version.
#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "matrix",
            fidelity = "matrix",
            tr_fun = "function"),
  function(absorption, fidelity, tr_fun, directions = 8, symm = TRUE, override, resistance) {

    if (!missing(resistance)) {
      warning("The resistance parameter is depcrecated. Use the data parameter instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data, absorption, fidelity, tr_fun = tr_fun, directions = directions, symm = symm, override = override))
  })

# TODO: stop-gap for parameter changes. Remove in future version.
#' @rdname samc
setMethod(
  "samc",
  signature(data = "missing",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            tr_fun = "function"),
  function(absorption, fidelity, tr_fun, directions = 8, symm = TRUE, override, resistance, latlon) {

    if (!missing(resistance)) {
      warning("The resistance parameter is depcrecated. Use the data parameter instead", call. = FALSE)
      data <- resistance
    } else {
      stop("Invalid arguments. The data parameter must be specified.", call. = FALSE)
    }

    return(samc(data, absorption, fidelity, tr_fun = tr_fun, directions = directions, symm = symm, latlon = latlon, override = override))
  })

            absorption = "missing",
            fidelity = "missing",
            tr_fun = "missing",
            p_mat = "dgCMatrix"),
  function(p_mat, override) {

    if (!missing(override))
      warning("The override parameter is deprecated. See the samc_opt() function instead.", call. = FALSE)

    r = nrow(p_mat)
    c = ncol(p_mat)

    if (c != r) stop("Matrix is not square", call. = FALSE)
    if (p_mat[r, c] != 1) stop("The last element must be 1", call. = FALSE)
    if (sum(p_mat[r,]) != 1) stop("Last row must be all zeros with a 1 in the last element", call. = FALSE)
    if (!isTRUE(all.equal(Matrix::rowSums(p_mat), rep(1, r), check.names = FALSE))) stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues

    if (is.null(rownames(p_mat))) rownames(p_mat) <- 1:r
    if (is.null(colnames(p_mat))) colnames(p_mat) <- 1:c

    rn <- rownames(p_mat)[-r]
    cn <- colnames(p_mat)[-r]

    if (!isTRUE(all.equal(rn, cn)))
      stop("The row and col names of the Q matrix must be identical", call. = FALSE)

    if (any(duplicated(rn)))
      stop("The row and col names of the Q matrix must be unique", call. = FALSE)


    print("Warning: Some checks for manually created P matrices are still missing:")
    print("1) Discontinuous data will not work with the cond_passage() function.")
    print("2) Every disconnected region of the graph must have at least one non-zero absorption value.")
    # TODO The clumps value is a placeholder and needs to be calculated as a safety check for the cond_passage() function
    samc_obj <- methods::new("samc",
                             p = p_mat,
                             source = "matrix",
                             map = raster::raster(matrix()),
                             clumps = 1,
                             override = FALSE)

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "missing",
            absorption = "missing",
            fidelity = "missing",
            tr_fun = "missing",
            p_mat = "matrix"),
  function(p_mat, override) {
    p <- as(p_mat, "dgCMatrix")

    return(samc(p_mat = p, override = override))
  })



# #' @rdname samc
# `samc_opt<-` <- function(x, option, value) {
#   if (!is(x, "samc")) stop("x must be an samc-class object created with the samc() function.", call. = FALSE)
#
#   if(option == "override") {
#     if (is.logical(value)) {
#       x@override <- value
#     } else {
#       stop("The override option must be set to TRUE/FALSE.")
#     }
#   } else {
#     stop("Invalid option specified.", call. = FALSE)
#   }
# }

