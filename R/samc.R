# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R check.R
NULL


#' Create an samc object
#'
#' Create an samc object that contains the absorbing Markov chain data
#'
#' This function is used to create a \code{\link{samc-class}} object. There are
#' multiple options for creating this object.
#'
#' \strong{Option 1: Raster or Matrix Maps}
#'
#' \emph{\code{samc(data = matrix, absorption = matrix, fidelity = matrix, model = list())}}
#'
#' \emph{\code{samc(data = SpatRaster, absorption = SpatRaster, fidelity = SpatRaster, model = list())}}
#'
#' \emph{\code{samc(data = RasterLayer, absorption = RasterLayer, fidelity = RasterLayer, model = list())}}
#'
#' The \code{\link{samc-class}} object can be created from a combination of
#' resistance (or conductance), absorption, and fidelity data. These different landscape data
#' inputs must be the same type (a matrix, SpatRaster, or RasterLayer), and have identical
#' properties, including dimensions, location of NA cells, and CRS (if using
#' raster inputs).
#'
#' The \code{data} and \code{absorption} inputs are always mandatory for this approach. The
#' \code{fidelity} input is optional. If the \code{fidelity} input is not provided, then it
#' is assumed that there is no site fidelity (i.e., individuals will always move
#' to an adjacent cell each time step).
#'
#' The \code{model} parameter is mandatory. It is used when calculating the values for
#' the transition matrix. \code{model} must be constructed as a list with a
#' transition function, the number of directions (4 or 8), and if the transition
#' function is symmetric (TRUE or FALSE; currently not used). Here is the template:
#' \code{list(fun = `function`, dir = `numeric`, sym = `logical`)}
#'
#' When using raster inputs, SpatRaster objects (from the terra package) are recommended
#' over RasterLayer objects (from the raster package). Internally, samc is using SpatRaster
#' objects, which means RasterLayer objects are being converted to SpatRaster objects,
#' which is a source of memory inefficiency.
#'
#' \strong{Option 2: P Matrix}
#'
#' \emph{\code{samc(data = matrix)}}
#'
#' \emph{\code{samc(data = dgCMatrix)}}
#'
#' The \code{data} parameter can be used alone to create a \code{\link{samc-class}} object
#' directly from a preconstructed P matrix. This matrix must be either a base R
#' matrix, or a sparse matrix (dgCMatrix format) from the Matrix package. It
#' must meet the following requirements:
#' \itemize{
#'   \item The number of rows must equal the number of columns (a square matrix)
#'   \item Total absorption must be a single column on the right-hand side of the matrix
#'   \item At the bottom of the matrix, there must be a row filled with 0's except
#'   for the last element (bottom-right of the matrix diagonal), which must be set to 1
#'   \item Every disconnected region of the matrix must have at least one non-zero
#'   absorbing value
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
#' \strong{Additional Information}
#'
#' Depending on the data used to construct the samc-class object, some metrics
#' may cause crashes. This is a result of the underlying P matrix having specific
#' properties that make some equations unsolvable. One known case is a P matrix
#' that represents a disconnected graph, which can lead to the \code{cond_passage()}
#' function crashing. In terms of raster/matrix inputs, a disconnected graph
#' occurs when one or more pixels/cells are unreachable from other pixels/cells
#' due to the presence of a full barrier made up of NA values. In a raster, these
#' may be obvious as islands but can be as inconspicuous as a rogue isolated
#' pixel. There may be other currently unknown situations that lead to unsolvable
#' metrics.
#'
#' Future work is planned towards identifying these issues during the creation of
#' the samc-class object and handling them appropriately to prevent inadvertent
#' crashes.
#'
#' \strong{Version 3 Changes}
#'
#' Support for creating samc-class objects from TransitionLayer objects was removed
#' so that the package is not dependent on gdistance.
#'
#' \strong{Version 2 Changes}
#'
#' Version 1.5.0 officially removed support for the deprecated \code{resistance}, \code{tr_fun},
#' \code{directions}, \code{p_mat}, \code{latlon}, and \code{override} arguments. Old
#' code will have to be updated to the new samc() function structure in order to work.
#'
#'
#' @param data A \code{\link[terra]{SpatRaster-class}} or \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}} or Matrix package dgCMatrix sparse matrix.
#' @param absorption A \code{\link[terra]{SpatRaster-class}} or \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param fidelity A \code{\link[terra]{SpatRaster-class}} or \code{\link[raster]{RasterLayer-class}} or \code{\link[base]{matrix}}
#' @param model A list with args for constructing a transition matrix.
#' @param options A list of options that changes how the samc behaves computationally.
#'
#' @return A \code{\link{samc-class}} object
#'
#' @example inst/examples/example.R
#'
#' @export

setGeneric(
  "samc",
  function(data, absorption, fidelity, model, options = NULL) {
    standardGeneric("samc")
  })


#' @rdname samc
setMethod(
  "samc",
  signature(data = "SpatRaster",
            absorption = "SpatRaster",
            fidelity = "SpatRaster",
            model = "list"),
  function(data, absorption, fidelity, model, options = NULL) {
    options = .validate_options(options)
    model = .validate_model(model, options$method)

    tr_fun <- model$fun
    directions <-model$dir
    symm <- model$sym

    # Make sure the input data all aligns
    check(c(data, fidelity, absorption))

    data_minmax = terra::minmax(data)
    if (data_minmax["min", 1] < 0) {
      stop("The data must not have values <= 0", call. = FALSE)
    }

    fid_minmax = terra::minmax(fidelity)
    if (fid_minmax["min", 1] < 0 || fid_minmax["max", 1] > 1) {
      stop("Fidelity values must be in range of 0-1", call. = FALSE)
    }

    fidabs_minmax = terra::minmax(fidelity + absorption)
    if (fidabs_minmax["max", 1] > 1) {
      stop("No cells can have fidelity + absoprtion > 1", call. = FALSE)
    }

    if (!(directions %in% c(4, 8))) {
      stop("directions must be set to either 4 or 8", call. = FALSE)
    }

    #abs_vec <- as.vector(absorption)
    #fid_vec <- as.vector(fidelity)

    abs_minmax = terra::minmax(absorption)

    if (abs_minmax["min", 1] < 0 || abs_minmax["max", 1] > 1) {
      stop("Absorption values must be in range of 0-1", call. = FALSE)
    } else if (abs_minmax["max", 1] == 0) {
      stop("At least one cell must have an absorption value > 0", call. = FALSE)
    }

    samc_obj <- methods::new("samc",
                             data = methods::new("samc_data",
                                                 f = new("dgCMatrix"),
                                                 t_abs = numeric(0)),
                             conv_cache = NULL,
                             model = model,
                             source = "SpatRaster",
                             map = data,
                             crw_map = NULL,
                             names = NULL,
                             clumps = -1,
                             override = options$override,
                             solver = options$method,
                             threads = options$threads,
                             .cache = new.env())

    # Check for "clumps"
    cl = terra::patches(samc_obj@map, directions = directions, zeroAsNA = TRUE, allowGaps = FALSE)
    samc_obj@clumps = sum(!is.na(terra::unique(cl)[, 1]))

    if (samc_obj@clumps > 1) {
      print("Warning: Input contains disconnected regions. This does not work with the cond_passage() metric.")

      temp_abs = absorption[[1]]
      temp_abs[temp_abs > 0] = 1
      temp_abs = temp_abs * cl

      if (!all(1:samc_obj@clumps %in% terra::unique(temp_abs)[, 1])) stop("All disconnected regions must have at least one non-zero absorption value", call. = FALSE)

      rm(temp_abs)
    }
    rm(cl)
    gc()
    # Create the transition matrix
    if (options$method %in% c("direct", "iter")) {
      if (model$name == "RW") {
        samc_obj@data@f = .rw(data, absorption, fidelity, tr_fun, directions, symm)
        gc()

        samc_obj@data@t_abs = as.vector(terra::values(absorption))[terra::cells(absorption)]
      } else if (model$name == "CRW") {
        warning("CRW support is currently experimental and may see input changes")

        if (terra::is.lonlat(data)) warning("CRW does not properly adjust turning angles for lonlat yet.")

        crw_list = .crw(data, absorption, fidelity, tr_fun, directions, symm, model)
        #assign("myvar", crw_list)
        samc_obj@data@f = crw_list$tr
        gc()

        samc_obj@data@t_abs = crw_list$abs
        samc_obj@crw_map = crw_list$crw

      } else {
        stop("Unexpected error involving model name. Please report with a minimum reproducible example.", call. = FALSE)
      }
    } else if (options$method == "conv") {
      if (unlist(terra::global(data, "isNA")) > 0) stop("Convolution method currently does not support rasters with NA data")

      if (terra::is.lonlat(data)) warning("Convolution method currently does not properly correct directions for lonlat")

      samc_obj@conv_cache = .convolution(data, absorption, fidelity, directions, symm, options$threads)
      samc_obj@data@t_abs = as.vector(terra::values(absorption))[terra::cells(absorption)]
    }


    # TODO Update to terra
    if(is.na(raster::projection(data)) && raster::xres(data) != raster::yres(data)) {
      warning("Raster cells are not square (number of columns/rows is not propotional to the spatial extents). There is no defined projection to account for this, so the geocorrection may lead to distortion if the intent was for the raster cells to represent a uniformly spaced grid.", call. = FALSE)
    }

    map_cells = terra::cells(samc_obj@map)
    samc_obj@map[map_cells] <- 1:length(map_cells)

    samc_obj@.cache$dgf = numeric(nrow(samc_obj@data@f))
    samc_obj@.cache$dgf_exists = FALSE

    samc_obj@.cache$sc = .solver_cache();

    return(samc_obj)
  })


#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "RasterLayer",
            model = "list"),
  function(data, absorption, fidelity, model, options = NULL) {

    data = rasterize(data)
    absorption = rasterize(absorption)
    fidelity = rasterize(fidelity)

    samc_obj = samc(data, absorption, fidelity, model = model, options = options)
    samc_obj@source = "RasterLayer"

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "SpatRaster",
            absorption = "SpatRaster",
            fidelity = "missing",
            model = "list"),
  function(data, absorption, model, options = NULL) {

    fidelity <- data
    fidelity[is.finite(fidelity)] <- 0

    return(samc(data, absorption, fidelity, model, options = options))
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "RasterLayer",
            absorption = "RasterLayer",
            fidelity = "missing",
            model = "list"),
  function(data, absorption, model,options = NULL) {

    data = rasterize(data)
    absorption = rasterize(absorption)

    samc_obj = samc(data, absorption, model = model, options = options)
    samc_obj@source = "RasterLayer"

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "matrix",
            model = "list"),
  function(data, absorption, fidelity, model, options = NULL) {

    data <- rasterize(data)
    absorption <- rasterize(absorption)
    fidelity <- rasterize(fidelity)

    samc_obj = samc(data, absorption, fidelity, model, options = options)
    samc_obj@source = "matrix"

    return(samc_obj)
  })



#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "matrix",
            fidelity = "missing",
            model = "list"),
  function(data, absorption, model, options = NULL) {

    data <- rasterize(data)
    absorption <- rasterize(absorption)

    samc_obj = samc(data, absorption, model = model, options = options)
    samc_obj@source = "matrix"

    return(samc_obj)
  })


#
# P matrix ----
#

#' @rdname samc
setMethod(
  "samc",
  signature(data = "dgCMatrix",
            absorption = "missing",
            fidelity = "missing",
            model = "missing"),
  function(data, options = NULL) {
    options = .validate_options(options)

    r = nrow(data)
    c = ncol(data)

    if (c != r) stop("Matrix is not square", call. = FALSE)
    if (data[r, c] != 1) stop("The last element must be 1", call. = FALSE)
    if (sum(data[r,]) != 1) stop("Last row must be all zeros with a 1 in the last element", call. = FALSE)
    if (!isTRUE(all.equal(Matrix::rowSums(data), rep(1, r), check.names = FALSE))) stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues


    q_mat <- methods::as(methods::as(data[-r, -c], "CsparseMatrix"), "generalMatrix")
    abs_total <- data[-r, c]

    if (!isTRUE(all.equal(Matrix::rowSums(q_mat) + abs_total, rep(1, length(abs_total)), check.names = FALSE))) stop("All row sums must be equal to 1", call. = FALSE) # Use all.equal() to avoid numerical precision issues

    if (is.null(rownames(q_mat)) & is.null(colnames(q_mat))) {
      nm = NULL
    } else if (!isTRUE(all.equal(rownames(q_mat), colnames(q_mat)))) {
      stop("The row and col names of the Q matrix must be identical", call. = FALSE)
    } else if (any(duplicated(rownames(q_mat)))) {
      stop("The row and col names of the Q matrix must be unique", call. = FALSE)
    } else {
      nm = rownames(q_mat)
      rownames(q_mat) = NULL
      colnames(q_mat) = NULL
    }

    q_mat@x = -q_mat@x
    Matrix::diag(q_mat) <- Matrix::diag(q_mat) + 1

    print("Warning: Some checks for manually created P matrices are still missing:")
    print("1) Discontinuous data will not work with the cond_passage() function.")
    print("2) Every disconnected region of the graph must have at least one non-zero absorption value.")
    # TODO The clumps value is a placeholder and needs to be calculated as a safety check for the cond_passage() function
    samc_obj <- methods::new("samc",
                             data = methods::new("samc_data",
                                                 f = q_mat,
                                                 t_abs = abs_total),
                             model = list(name = "RW"),
                             crw_map = NULL,
                             source = "transition",
                             map = terra::rast(),
                             names = nm,
                             clumps = -1,
                             threads = options$threads,
                             override = options$override,
                             solver = options$method)

    samc_obj@.cache$dgf = numeric(nrow(samc_obj@data@f))
    samc_obj@.cache$dgf_exists = FALSE

    samc_obj@.cache$sc = .solver_cache();

    return(samc_obj)
  })

#' @rdname samc
setMethod(
  "samc",
  signature(data = "matrix",
            absorption = "missing",
            fidelity = "missing",
            model = "missing"),
  function(data, options = NULL) {
    p <- methods::as(methods::as(data, "CsparseMatrix"), "generalMatrix")

    return(samc(data = p, options = options))
  })
