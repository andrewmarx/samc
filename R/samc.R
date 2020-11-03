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
#' The resistance and absorption inputs are always mandatory, whereas the
#' fidelity input is optional. If the fidelity input is not provided, then it it
#' is assumed that there is no site fidelity (i.e., individuals will always move
#' to an adjacent cell each time step).
#'
#' The latlon parameter is required if the landscape data inputs are RasterLayer
#' objects. The package does not attempt to determine this automatically, and it
#' does not assume a default. Users must set it to TRUE if they are using
#' latitude and longitude data.
#'
#' The tr_fun parameter is mandatory. It used when calculating the values for
#' the transition matrix. Internally, this is passed to the \code{\link[gdistance]{transition}}
#' function in the gdistance package to create the transition matrix.
#'
#' \strong{Option 2: P Matrix Input}
#'
#' The p_mat parameter can be used to create a \code{\link{samc-class}} object
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
#' \strong{Other Parameters}
#'
#' The override parameter is optional. To prevent users from unintentionally
#' running memory intensive versions of functions that could make their systems
#' non-responsive or crash software, it is set to FALSE by default. For various
#' reasons, it can be set to TRUE. In particular, a user might do this if they
#' are using a very small landscape dataset, or perhaps for a moderately sized
#' dataset if they have access to a system with exceptionally large amounts of
#' RAM. Before setting this to TRUE, users should read the Performance vignette/
#' article to understand the expected memory requirements. They should also
#' consider starting with scaled down version of their data and then gradually
#' scaling back up while monitoring their memory usage as a means to gauge what
#' is reasonable for their system.
#'
#' \strong{Additional Information}
#'
#' Depending on the data used to construct the samc-class object, some metrics
#' may cause crashes. This is a result of the underlying P matrix having specific
#' properties that make some equations unsolvable. One known case is a P matrix
#' that represents a disconnected graph, which can lead to the cond_passage()
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
#' @param resistance A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param absorption A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param fidelity A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param latlon Logical (\code{TRUE} or \code{FALSE}) indicating whether the rasters use latitude/longitude
#' @param tr_fun A function to calculate the transition values in the \code{\link[gdistance]{transition}} function
#' @param p_mat A base R \code{\link[base]{matrix}} object or Matrix package dgCMatrix sparse matrix
#' @param override Optional flag to prevent accidentally running memory intensive functions. Defaults to \code{FALSE}
#' @param ... Placeholder
#'
#' @return A spatial absorbing Markov chain object
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "samc",
  function(resistance, absorption, fidelity, latlon, tr_fun, p_mat, ...) {
    standardGeneric("samc")
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            latlon = "logical",
            tr_fun = "function",
            p_mat = "missing"),
  function(resistance, absorption, fidelity, latlon, tr_fun, override = FALSE) {

    if (!is.logical(override))
      stop("The override parameter must be set to TRUE or FALSE")

    # Make sure the input data all aligns
    check(resistance, absorption)
    check(resistance, fidelity)

    if (any(resistance[] <= 0, na.rm = TRUE)) {
      stop("The resistance data must not have values <= 0")
    }

    if (any(absorption[] < 0, na.rm = TRUE)) {
      stop("The absorption data must not have values <= 0")
    }

    if (any(absorption[] > 1, na.rm = TRUE)) {
      stop("The absorption data must not have values > 1")
    }

    if (sum(absorption[], na.rm = TRUE) == 0) {
      stop("At least one cell must have an absorption value > 0")
    }

    if (any(fidelity[] < 0, na.rm = TRUE)) {
      stop("The fidelity data must not have values < 0")
    }

    if (any(fidelity[] > 1, na.rm = TRUE)) {
      stop("The fidelity data must not have values > 1")
    }

    if (any((fidelity[] + absorption[]) > 1, na.rm = TRUE)) {
      stop("No cells can have fidelity + absoprtion > 1")
    }


    # Create map template
    m <- resistance
    m[] <- is.finite(m[])

    # Check for "clumps"
    cl <- raster::clump(m, directions = 8, gaps = FALSE)
    clumps <- sum(!is.na(unique(cl[])))

    if (clumps > 1) {
      print("Warning: Input contains disconnected regions. This does not work with the cond_passage() metric.")

      temp_abs <- absorption
      temp_abs[temp_abs > 0] <- 1
      temp_abs <- temp_abs * cl

      if (!all(1:clumps %in% unique(temp_abs[]))) stop("All disconnected regions must have at least one non-zero absorption value")
    }


    # Create the transition matrix
    tr <- gdistance::transition(resistance, transitionFunction = tr_fun, 8)
    if (latlon) {
      tr <- gdistance::geoCorrection(tr, type = "c")
    } else {
      tr <- gdistance::geoCorrection(tr)
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

    # Assemble final

    samc_mat <- methods::new("samc", p = p, source = "map", map = m, clumps = clumps, override = override)

    return(samc_mat)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            latlon = "logical",
            tr_fun = "function",
            p_mat = "missing"),
  function(resistance, absorption, latlon, tr_fun, override = FALSE) {

    fidelity <- resistance
    fidelity[is.finite(fidelity)] <- 0

    return(samc(resistance, absorption, fidelity, latlon, tr_fun, override = override))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "matrix",
            absorption = "matrix",
            fidelity = "matrix",
            latlon = "missing",
            tr_fun = "function",
            p_mat = "missing"),
  function(resistance, absorption, fidelity, tr_fun, override = FALSE) {

    resistance <- raster::raster(resistance, xmn = 0.5, xmx = ncol(resistance) + 0.5, ymn = 0.5, ymx = nrow(resistance) + 0.5)
    absorption <- raster::raster(absorption, xmn = 0.5, xmx = ncol(absorption) + 0.5, ymn = 0.5, ymx = nrow(absorption) + 0.5)
    fidelity <- raster::raster(fidelity, xmn = 0.5, xmx = ncol(fidelity) + 0.5, ymn = 0.5, ymx = nrow(fidelity) + 0.5)

    #fidelity[is.finite(fidelity)] <- 0

    return(samc(resistance, absorption, fidelity, FALSE, tr_fun, override = override))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "matrix",
            absorption = "matrix",
            fidelity = "missing",
            latlon = "missing",
            tr_fun = "function",
            p_mat = "missing"),
  function(resistance, absorption, tr_fun, override = FALSE) {

    resistance <- raster::raster(resistance, xmn = 0.5, xmx = ncol(resistance) + 0.5, ymn = 0.5, ymx = nrow(resistance) + 0.5)
    absorption <- raster::raster(absorption, xmn = 0.5, xmx = ncol(absorption) + 0.5, ymn = 0.5, ymx = nrow(absorption) + 0.5)

    return(samc(resistance, absorption, latlon = FALSE, tr_fun = tr_fun, override = override))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "missing",
            absorption = "missing",
            fidelity = "missing",
            latlon = "missing",
            tr_fun = "missing",
            p_mat = "dgCMatrix"),
  function(p_mat, override = FALSE) {
    r = nrow(p_mat)
    c = ncol(p_mat)

    if (c != r) stop("Matrix is not square")
    if (p_mat[r, c] != 1) stop("The last element must be 1")
    if (sum(p_mat[r,]) != 1) stop("Last row must be all zeros with a 1 in the last element")
    if (!isTRUE(all.equal(Matrix::rowSums(p_mat), rep(1, r)))) stop("All row sums must be equal to 1") # Use all.equal() to avoid numerical precision issues

    print("Warning: Some checks for manually created P matrices are still missing:")
    print("1) Discontinuous data will not work with the cond_passage() function.")
    print("2) Every disconnected region of the graph must have at least one non-zero absorption value.")
    # TODO The clumps value is a placeholder and needs to be calculated as a safety check for the cond_passage() function
    samc_obj <- methods::new("samc",
                             p = p_mat,
                             source = "matrix",
                             map = raster::raster(matrix()),
                             clumps = 1,
                             override = override)

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(resistance = "missing",
            absorption = "missing",
            fidelity = "missing",
            latlon = "missing",
            tr_fun = "missing",
            p_mat = "matrix"),
  function(p_mat, override = FALSE) {
    p <- as(p_mat, "dgCMatrix")

    return(samc(p_mat = p, override = override))
  })
