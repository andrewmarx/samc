# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.


#' data class
#'
#' Contains the data fields used in the samc-class
#'
#' @noRd

setClass(
  # set the name of the class
  "samc_data",

  # define the slots
  slots = list(q = "dgCMatrix",
               r = "matrix")

  # set default values
  #prototype = list(p = NA)

  # create a function to validate the data
  # validity=function(object)
  # {
  #   return(TRUE)
  # }
)
