# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R check.R
NULL


#' Create an samc object
#'
#' Create an samc object that contains the absorbing Markov chain data
#'
#' This function is used to create a \code{\link{samc-class}} object from
#' landscape data. Some of the inputs are mandatory, whereas others are
#' optional. The different landscape data inputs must be the same type (a matrix
#' or RasterLayer), and have identical properties, including dimensions,
#' location of NA cells, and CRS (if using RasterLayers).
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
#'
#' @param resistance A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param absorption A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param fidelity A \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param latlon Logical (\code{TRUE} or \code{FALSE}) indicating whether the rasters use latitude/longitude
#' @param tr_fun A function to calculate the transition values in the \code{\link[gdistance]{transition}} function
#' @param p_mat An option to provide the P matrix directly
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
    fid_norm = FALSE

    if (!is.logical(override))
      stop("The override parameter must be set to TRUE or FALSE")

    # Make sure the input data all aligns
    check(resistance, absorption)
    check(resistance, fidelity)

    if (any(resistance[] <= 0, na.rm = TRUE)) {
      stop("The resistance data must not have values <= 0")
    }

    if (any(absorption[] <= 0, na.rm = TRUE)) {
      stop("The absorption data must not have values <= 0")
    }

    if (any(absorption[] > 1, na.rm = TRUE)) {
      stop("The absorption data must not have values > 1")
    }

    if (any(fidelity[] < 0, na.rm = TRUE)) {
      stop("The fidelity data must not have values < 0")
    }

    if (any(fidelity[] > 1, na.rm = TRUE)) {
      stop("The fidelity data must not have values > 1")
    }

    tr <- gdistance::transition(resistance, transitionFunction = tr_fun, 8)
    if (latlon) {
      tr <- gdistance::geoCorrection(tr, type = "c")
    } else {
      tr <- gdistance::geoCorrection(tr)
    }

    tr_mat <- gdistance::transitionMatrix(tr)

    abs_vec <- as.vector(absorption)
    fid_vec <- as.vector(fidelity)
    tr_mat <- methods::as(tr_mat, "dgTMatrix")

    # 'Remove' the null values by changing the index values using a lookup vector
    # Originally created later in the code, but addition of the fid_norm code
    # interfered by causing the dgTMatrix to populate i values
    lookup_vec <- 0:(length(unique(tr_mat@i)))
    names(lookup_vec) <- c(sort(unique(tr_mat@i)), length(abs_vec))

    if (fid_norm) { # TODO: Needs to be finished. Note that fid_norm is set to false above
      diag(tr_mat) <- fid_vec
      tr_mat@x <- (1 - abs_vec[tr_mat@i + 1]) * tr_mat@x / Matrix::rowSums(tr_mat)[tr_mat@i + 1]
    } else {
      diag(tr_mat) <- 0
      tr_mat@x <- (1 - abs_vec[tr_mat@i + 1] - fid_vec[tr_mat@i + 1]) * tr_mat@x / Matrix::rowSums(tr_mat)[tr_mat@i + 1]
      diag(tr_mat) <- fid_vec
    }


    # Combine the transition matrix with the absorbing data
    samc_df <- data.frame(i = c(tr_mat@i, (0:(length(abs_vec) - 1))[is.finite(abs_vec)], length(abs_vec)),
                         j = c(tr_mat@j, rep(length(abs_vec), sum(is.finite(abs_vec))), length(abs_vec)),
                         x = c(tr_mat@x, abs_vec[is.finite(abs_vec)], 1))

    samc_df$i <- lookup_vec[as.character(samc_df$i)]
    samc_df$j <- lookup_vec[as.character(samc_df$j)]

    samc_df <- samc_df[!is.na(samc_df$i),]
    samc_df <- samc_df[!is.na(samc_df$j),]

    m <- resistance
    m[] <- is.finite(m[])

    # Create the final sparse matrix
    samc_mat <- methods::new("samc",
                    p = Matrix::sparseMatrix(i = samc_df$i,
                                             j = samc_df$j,
                                             x = samc_df$x,
                                             index1 = FALSE),
                    source = "map",
                    map = m,
                    override = override)

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

    resistance <- raster::raster(resistance, xmn = 1, xmx = ncol(resistance), ymn = 1, ymx = nrow(resistance))
    absorption <- raster::raster(absorption, xmn = 1, xmx = ncol(absorption), ymn = 1, ymx = nrow(absorption))
    fidelity <- raster::raster(fidelity, xmn = 1, xmx = ncol(fidelity), ymn = 1, ymx = nrow(fidelity))

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

    resistance <- raster::raster(resistance, xmn = 1, xmx = ncol(resistance), ymn = 1, ymx = nrow(resistance))
    absorption <- raster::raster(absorption, xmn = 1, xmx = ncol(absorption), ymn = 1, ymx = nrow(absorption))

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

    samc_obj <- methods::new("samc",
                             p = p_mat,
                             source = "matrix",
                             map = raster::raster(matrix()),
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
