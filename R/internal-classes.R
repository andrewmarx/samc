# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.
NULL

#' data class
#'
#' Contains the data fields used in the samc-class
#'
#' @slot q Q matrix
#' @slot r R matrix
#' @slot t_abs Total absorption
#'
#' @name samc_data-class
#' @keywords internal

setClass(
  # set the name of the class
  "samc_data",

  # define the slots
  slots = list(q = "dgCMatrix",
               r = "matrix",
               t_abs = "numeric")

  # set default values
  #prototype = list(p = NA)

  # create a function to validate the data
  # validity=function(object)
  # {
  #   return(TRUE)
  # }
)
