# Copyright (c) 2020 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.
#
# This file is for internal functions. They are subject to change and should not
# be used by users.
#

#' Validate time steps
#'
#' Performs several checks to make sure a vector of time steps is valid
#'
#' @param x A vector object to be validated as time steps
#' @noRd
.validate_time_steps <- function(x) {
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

  if (any(x > 10000))
    stop("Due to how the short-term metrics are calculated and the way that
  decimal numbers are handled by computers, numerical issues related to
  precision arise when a time step value is too high. Currently, a hard
  limit of 10000 time steps is enforced to encourage users to more
  seriously consider how many time steps are relevant to their use case.
  For example, if a single time step represents 1 day, then the current
  limit represents 24.7 years. There is flexibility to increase the limit
  if a justification can be made for it, but it's far more likely that
  users will generally want far fewer time steps for ecologically relevant
  results and to avoid the cummulative precision issues.", call. = FALSE)
}


#' Validate location vectors
#'
#' Performs several checks to make sure a vector locations is valid
#'
#' @param samc samc-class object
#' @param x A vector object to be validated as locations
#' @noRd
.validate_locations <- function(samc, x) {
  if (!is.numeric(x))
    stop("Locations must be a positive integer or a vector of positive integers", call. = FALSE)

  if (sum(is.na(x)) > 0)
    stop("NA values are not valid locations", call. = FALSE)

  if (any(x %% 1 != 0))
    stop("Decimal values are not valid locations", call. = FALSE)

  if (any(x < 1))
    stop("All location values must be positive (greater than 0)", call. = FALSE)

  if (any(x > (nrow(samc@p) - 1)))
    stop("Location values cannot exceed the number of nodes in the landscape", call. = FALSE)
}


#' Rasterize matrices
#'
#' Convert a matrix to a RasterLayer. Ensures consistency of conversion throughout the package
#'
#' @param x A matrix
#' @noRd
.rasterize <- function(x) {
  return(raster::raster(x, xmn = 0.5, xmx = ncol(x) + 0.5, ymn = 0.5, ymx = nrow(x) + 0.5))
}
