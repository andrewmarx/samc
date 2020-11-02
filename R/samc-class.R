# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.


#' samc class
#'
#' S4 class to manage SAMC data.
#'
#' The samc class is used to help ensure that the package is used correctly and
#' to minimize the possibility for users to accidentally produce nonsensical
#' results that may not be obviously incorrect. This class contains the p matrix
#' necessary for all the calculations in the package, and enforces its type so
#' that users are less likely to inadvertently alter it in a way that will cause
#' issues in calculations.
#'
#' The class also contains a RasterLayer object derived from the input data.
#' This object is used for checking inputs and mapping vector data in other
#' functions.
#'
#' Finally, an override flag is used to help ensure that users do not
#' accidentally run memory intensive versions of functions that can cause their
#' systems to become non-responsive or for software to crash.
#'
#' The \code{\link{samc}} function is used to create \code{\link{samc-class}}
#' objects.
#'
#' @slot p The transition probability matrix \emph{P}.
#' @slot source Information about the data source for the P matrix
#' @slot map Used to verify landscape inputs and mapping of vector data.
#' @slot clumps Number of discontinuous regions in data
#' @slot override Used to prevent accidental use of memory intensive functions.


setClass(
  # set the name of the class
  "samc",

  # define the slots
  slots = list(p = "dgCMatrix",
               source = "character",
               map = "RasterLayer",
               clumps = "numeric",
               override = "logical")

  # set default values
  #prototype = list(p = NA)

  # create a function to validate the data
  # validity=function(object)
  # {
  #   return(TRUE)
  # }
  )
