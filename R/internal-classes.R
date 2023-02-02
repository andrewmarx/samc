# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.
NULL

#' data class
#'
#' Contains the data fields used in the samc-class
#'
#' @slot f F matrix
#' @slot t_abs Total absorption
#' @slot c_abs Component absorption states
#'
#' @name samc_data-class
#' @keywords internal

setClass(
  # set the name of the class
  "samc_data",

  # define the slots
  slots = list(f = "CsparseMatrix",
               t_abs = "numeric",
               c_abs = "matrix")

  # set default values
  #prototype = list(p = NA)

  # create a function to validate the data
  # validity=function(object)
  # {
  #   return(TRUE)
  # }
)


#' samc char_null class
#'
#' Class for grouping character and NULL data types
#'
#' @name char_null-class
#' @keywords internal
#'
setClassUnion("char_null", c("character", "NULL"))
