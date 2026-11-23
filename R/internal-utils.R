# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

# This file is for internal functions. They are subject to change and should not
# be used by users.


#' Test whether a samc object uses a random walk (RW) state space
#'
#' TRUE for models with one transient state per cell. SSF models only differ
#' from RW in how the transition values are calculated, so they are included.
#'
#' @param samc samc model
#' @noRd
.is_rw = function(samc) {
    samc@model$name %in% c("RW", "SSF")
}

#' Test whether a samc object uses the correlated random walk (CRW) model
#'
#' @param samc samc model
#' @noRd
.is_crw = function(samc) {
    samc@model$name == "CRW"
}

#' Resolve per-state total absorption, collapsing CRW substates
#'
#' Returns the total absorption vector. For CRW models the per-substate values
#' are summarized to one value per state using the mean.
#'
#' @param samc samc model
#' @noRd
.resolve_t_abs = function(samc) {
    if (.is_crw(samc)) .summarize_crw(samc, samc@data@t_abs, mean) else samc@data@t_abs
}


#' Used to disable CRW
#'
#' Disable CRW
#'
#' @param samc samc model
#' @noRd
.disable_crw = function(samc) {
    if (.is_crw(samc)) stop("Metric/parameter combination not currently supported for CRW", call. = FALSE)
}

#' Used to disable convolution
#'
#' Disable convolution
#'
#' @param samc samc model
#' @noRd
.disable_conv = function(samc) {
    if (samc@solver == "conv") stop("Metric/parameter combinaton not currently supported for the convolution algorithm", call. = FALSE)
}
