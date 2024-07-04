# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.


#' location class
#'
#' Union class for location inputs
#'
#' The location class is a union class of the "numeric" and "character" classes.
#' Users generally do not need to worry about it except to know that any method
#' parameter with "location" as the type can have either an integer or a character
#' name provided as input.
#'
#' @export location

setClassUnion("location", c("numeric", "matrix", "character"))
