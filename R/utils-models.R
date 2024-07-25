# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#' Validate model
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_model <- function(x, method) {
  names <- names(x)

  dup_args <- names[duplicated(names)]
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
.validate_rw <- function(x, method) {
  names <- names(x)

  args = c("name", "fun", "dir", "sym")
  methods = c("direct", "iter", "conv")

  missing_args <- args[!(args %in% names)]
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

  unknown_args <- names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)
}

#' Validate transition args for RW
#'
#' Validates the model for the samc() function
#'
#' @param x A list
#' @noRd
.validate_crw <- function(x, method) {
  names <- names(x)

  args = c("name", "fun", "dir", "sym", "dist", "kappa")
  methods = c("direct", "iter")


  missing_args <- args[!(args %in% names)]
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

  unknown_args <- names[!(names %in% args)]
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
.validate_ssf <- function(x, method) {
  names <- names(x)

  args = c("name", "fun", "dir", "sym", "ssc")
  methods = c("direct", "iter")

  missing_args <- args[!(args %in% names)]
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

  unknown_args <- names[!(names %in% args)]
  if (length(unknown_args) > 0)
    stop(paste("Unknown argument in model:", unknown_args), call. = FALSE)

  if (!is(x$ssc, "numeric"))
    stop("ssc must be single numeric value.", call. = FALSE)

  if (length(x$ssc) != 1)
    stop("ssc must be single numeric value.", call. = FALSE)

  if (!is.finite(x$ssc))
    stop("ssc must be single numeric value.", call. = FALSE)
}


