# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.

#' Is single pos int
#'
#' Validates that a variable is a single positive integer
#'
#' @param x anything
#' @noRd
.is_single_positive_integer = function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && x == as.integer(x) && x > 0
}


#' Is single string
#'
#' Validates that a variable is a single non-NA character string
#'
#' @param x anything
#' @noRd
.is_single_string = function(x) {
    is.character(x) && length(x) == 1 && !is.na(x)
}


#' Is single logical
#'
#' Validates that a variable is a single non-NA logical
#'
#' @param x anything
#' @noRd
.is_single_logical = function(x) {
    is.logical(x) && length(x) == 1 && !is.na(x)
}


#' Is valid dir
#'
#' Validates that a variable is a single number, without dimensions, equal to
#' either 4 or 8
#'
#' @param x anything
#' @noRd
.is_valid_dir = function(x) {
    is.numeric(x) && is.null(dim(x)) && length(x) == 1 && !is.na(x) && x %in% c(4, 8)
}


#' Validate model
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_model = function(x, method) {
    names = names(x)

    if (length(x) > 0 && (is.null(names) || any(is.na(names) | names == "")))
        stop("model: all elements must be named", call. = FALSE)

    dup_args = unique(names[duplicated(names)])
    if (length(dup_args) > 0)
        stop(paste("Duplicate argument in model:", paste(dup_args, collapse = ", ")), call. = FALSE)

    if (!("name" %in% names)) {
        x$name = "RW"
    }

    switch(
        as.character(x$name),
        RW = .validate_rw(x, method),
        CRW = .validate_crw(x, method),
        SSF = .validate_ssf(x, method),
        stop("Invalid model name", call. = FALSE)
    )

    return(x)
}


#' Validate args common to all transition models
#'
#' Checks the shared argument structure (required args, fun/dir/sym types,
#' method, and unknown args) used by every model validator.
#'
#' @param x A list
#' @param args Character vector of allowed argument names
#' @param funs Character vector of the named functions supported by the model
#' @param methods Character vector of allowed methods
#' @param method The method to validate against
#' @noRd
.validate_model_common = function(x, args, funs, methods, method) {
    names = names(x)

    missing_args = args[!(args %in% names)]
    if (length(missing_args) > 0)
        stop(paste("Missing argument in model:", paste(missing_args, collapse = ", ")), call. = FALSE)

    if (!is.function(x$fun)) {
        if (!.is_single_string(x$fun))
            stop("`fun` must be a function or the name of a supported named function", call. = FALSE)

        if (!(x$fun %in% funs)) {
            supported = paste0("\"", funs, "\"", collapse = ", ")
            stop(paste0("`fun` = \"", x$fun, "\" is not a named function supported by the ", x$name,
                        " model. Supported names: ", supported), call. = FALSE)
        }
    }

    if (!.is_valid_dir(x$dir))
        stop("`dir` must be set to either 4 or 8", call. = FALSE)

    if (!.is_single_logical(x$sym))
        stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)

    if (!(method %in% methods))
        stop("Invalid method for model", call. = FALSE)

    unknown_args = names[!(names %in% args)]
    if (length(unknown_args) > 0)
        stop(paste("Unknown argument in model:", paste(unknown_args, collapse = ", ")), call. = FALSE)
}


#' Validate transition args for RW
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_rw = function(x, method) {
    args = c("name", "fun", "dir", "sym")
    funs = c("1/mean(x)")
    methods = c("direct", "iter", "conv")

    .validate_model_common(x, args, funs, methods, method)

    if (method == "conv" && !identical(x$fun, "1/mean(x)"))
        stop("Convolution currently only supports the '1/mean(x)' named function.", call. = FALSE)
}


#' Validate transition args for CRW
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_crw = function(x, method) {
    args = c("name", "fun", "dir", "sym", "dist", "kappa")
    funs = c("1/mean(x)")
    methods = c("direct", "iter")

    .validate_model_common(x, args, funs, methods, method)

    if (x$dist == "vonMises") {
        if (is(x$kappa, "SpatRaster")) {
            if (!terra::is.num(x$kappa)) {
                stop("raster kappa must be numeric.")
            }

            k_minmax = terra::minmax(x$kappa)
            if (k_minmax["min", 1] < 0) {
                stop("raster kappa cannot have negative values.", call. = FALSE)
            }

            if (k_minmax["max", 1] == Inf) {
                stop("raster kappa cannot have Inf values.", call. = FALSE)
            }

        } else if (is.numeric(x$kappa)) {
            if (length(x$kappa) != 1)
                stop("scalar kappa must be single non-negative numeric value.", call. = FALSE)

            if (!is.finite(x$kappa))
                stop("scalar kappa must be single non-negative numeric value.", call. = FALSE)

            if (x$kappa < 0)
                stop("scalar kappa must be single non-negative numeric value.", call. = FALSE)
        } else {
            stop("kappa must be a non-negative scalar or raster map")
        }
    } else {
        stop(paste("Invalid distribution name:", x$dist), call. = FALSE)
    }
}


#' Validate transition args for SSF
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_ssf = function(x, method) {
    args = c("name", "fun", "dir", "sym", "ssc")
    funs = c("x[2]")
    methods = c("direct", "iter")

    .validate_model_common(x, args, funs, methods, method)

    if (!is(x$ssc, "numeric"))
        stop("ssc must be single numeric value.", call. = FALSE)

    if (length(x$ssc) != 1)
        stop("ssc must be single numeric value.", call. = FALSE)

    if (!is.finite(x$ssc))
        stop("ssc must be single numeric value.", call. = FALSE)
}


#' Validate time steps
#'
#' Performs several checks to make sure a vector of time steps is valid
#'
#' @param x A vector object to be validated as time steps
#' @noRd
.validate_time_steps = function(x) {
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
    #     stop("Due to how the short-term metrics are calculated and the way that
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
.validate_locations = function(samc, x) {
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
.validate_names = function(vec, x) {
    invalid_names = x[!(x %in% vec)]

    if (length(invalid_names > 0)) {
        #print(vec)
        #print(x)
        stop(paste("\nInvalid location name:", invalid_names), call. = FALSE)
    }
}


#' Validate options
#'
#' Validates the options args for the samc() function
#'
#' @param x A list
#' @noRd
.validate_options = function(x) {
    opt_names = c("method", "threads", "override", "precision")

    if (is.null(x)) {
        x = list()
    }

    if (is.list(x)) {
        if (is.null(x$method)) x$method = "direct"
        if (is.null(x$threads)) x$threads = 1
        if (is.null(x$override)) x$override = FALSE
        if (is.null(x$precision)) x$precision = "double"
    } else {
        stop("options argument must be a list or left empty for default values", call. = FALSE)
    }

    x_names = names(x)
    if (any(is.na(x_names) | x_names == "")) {
        stop("options: all elements must be named", call. = FALSE)
    }

    inv_names = x_names[!(x_names %in% opt_names)]
    if (length(inv_names) > 0) {
        stop(paste("Invalid option names:", paste(inv_names, collapse = " ")), call. = FALSE)
    }

    if (!(.is_single_string(x$method) && x$method %in% c("direct", "iter", "conv"))) {
        stop("options: method must be 'direct', 'iter', or 'conv'", call. = FALSE)
    }
    if (!.is_single_positive_integer(x$threads)) stop("options: threads must be a single positive integer", call. = FALSE)
    if (!.is_single_logical(x$override)) stop("options: override must be TRUE or FALSE", call. = FALSE)
    if (!(.is_single_string(x$precision) && x$precision %in% c("single", "double"))) {
        stop("options: precision must be 'single' or 'double'", call. = FALSE)
    }

    # Single precision is only implemented for the convolution method
    if (x$precision == "single" && x$method != "conv") {
        stop(paste0("method: '", x$method, "' only supports precision: 'double'"), call. = FALSE)
    }

    x
}
