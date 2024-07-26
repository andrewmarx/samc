# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.


#' Validate model
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_model = function(x, method) {
  names = names(x)

  dup_args = names[duplicated(names)]
  if (length(dup_args) > 0)
    stop(paste("Duplicate argument in model:", dup_args), call. = FALSE)

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


#' Validate transition args for CRW
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_rw = function(x, method) {
  names = names(x)

  args = c("name", "fun", "dir", "sym")
  methods = c("direct", "iter", "conv")

  missing_args = args[!(args %in% names)]
  if (length(missing_args) > 0)
    stop(paste("Missing argument in model:", missing_args), call. = FALSE)

  if (!(is(x$fun, "function") || is(x$fun, "character"))) {
    stop("'fun' must be a supported named function or a user defined function")
  } else if (!(x$dir %in% c(4, 8))) {
    stop("`dir` must be set to either 4 or 8", call. = FALSE)
  } else if (!is(x$sym, "logical")) {
    stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)
  }

  if (!(method %in% methods))
    stop("Invalid method for model", call. = FALSE)

  if (method == "conv") {
    if (!is(x$fun, "character")) {
      stop("Convolution currently only supports the '1/mean(x)' named function.", call. = FALSE)
    } else if (x$fun != "1/mean(x)") {
      stop("Convolution currently only supports the '1/mean(x)' named function.", call. = FALSE)
    }
  }

  unknown_args = names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)
}


#' Validate transition args for RW
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_crw = function(x, method) {
  names = names(x)

  args = c("name", "fun", "dir", "sym", "dist", "kappa")
  methods = c("direct", "iter")

  missing_args = args[!(args %in% names)]
  if (length(missing_args) > 0)
    stop(paste("Missing argument in model:", missing_args), call. = FALSE)

  if (!(is(x$fun, "function") || is(x$fun, "character"))) {
    stop("'fun' must be a supported named function or a user defined function")
  } else if (!(x$dir %in% c(4, 8))) {
    stop("`dir` must be set to either 4 or 8", call. = FALSE)
  } else if (!is(x$sym, "logical")) {
    stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)
  }

  if (!(method %in% methods))
    stop("Invalid method for model", call. = FALSE)

  unknown_args = names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)


  if (x$dist == "vonMises") {
    if (!is(x$kappa, "numeric"))
      stop("kappa must be single non-negative numeric value.", call. = FALSE)

    if (length(x$kappa) != 1)
      stop("kappa must be single non-negative numeric value.", call. = FALSE)

    if (!is.finite(x$kappa))
      stop("kappa must be single non-negative numeric value.", call. = FALSE)

    if (x$kappa < 0)
      stop("kappa must be single non-negative numeric value.", call. = FALSE)
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
  names = names(x)

  args = c("name", "fun", "dir", "sym", "ssc")
  methods = c("direct", "iter")

  missing_args = args[!(args %in% names)]
  if (length(missing_args) > 0)
    stop(paste("Missing argument in model:", missing_args), call. = FALSE)

  if (!(is(x$fun, "function") || is(x$fun, "character"))) {
    stop("'fun' must be a supported named function or a user defined function")
  } else if (!(x$dir %in% c(4, 8))) {
    stop("`dir` must be set to either 4 or 8", call. = FALSE)
  } else if (!is(x$sym, "logical")) {
    stop("`sym` must be set to either TRUE or FALSE", call. = FALSE)
  }

  if (!(method %in% methods))
    stop("Invalid method for model", call. = FALSE)

  unknown_args = names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)

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

  if (length(invalid_names > 0)){
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
