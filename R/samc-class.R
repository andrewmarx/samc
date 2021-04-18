# Copyright (c) 2019 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include internal-classes.R
NULL

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
#' The \code{\link{samc}} function is used to create \code{\link{samc-class}}
#' objects.
#'
#' The samc-class slots are subject to change, so users should not be using the
#' \code{@} operator to access or change them. Doing so leads to the risk of broken
#' code in the future. Instead, where relevant, the \code{$} operator can be used
#' to get and set components of the class safely. This is a current list of what can
#' be accessed and modified in the class:
#'
#' \itemize{
#'   \item \strong{override}
#'
#'   Some analyses are memory intensive and have the potential to make a user's
#'   system non-responsive or crash. By default, a samc-class object cannot be used
#'   in these analyses to prevent unintentional loss of work. In some cases, users
#'   may wish to use these particular analyses, in which case this behavior can
#'   be overridden. To get the current state of the override, use \code{samc_obj$override}.
#'   To enable the use of the analyses, the override can be set to \code{TRUE} using
#'   \code{samc_obj$override <- TRUE}. Before enabling the override, users should
#'   familiarize themselves with the Performance vignette. They should also consider
#'   starting with scaled down versions of their data and then gradually scale
#'   back up while monitoring their memory usage to guage what their system can
#'   handle.
#'
#'   \item \strong{q_matrix}
#'
#'   Advanced users may wish to have direct access to the Q matrix for developing
#'   custom calculations/analyses. Assumptions should not be made about the internal
#'   structure of the P and Q matrices in the samc-class, nor should it be assumed
#'   that they will not change in the future. To safely access the Q matrix, use
#'   \code{samc_obj$q_matrix}. The Q matrix inside of the samc-class cannot be
#'   modified.
#' }
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
