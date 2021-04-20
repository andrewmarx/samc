# Copyright (c) 2021 Andrew Marx. All rights reserved.
# Licensed under GPLv3.0. See LICENSE file in the project root for details.

#' @include samc-class.R
NULL

#' Access samc-class components
#'
#' Allows users to access a subset of the samc-class components
#'
#' @param x samc-class object
#' @param name Component of the samc-class to access
#'
#' @return Component value
#'
#' @name samc-class-access
#' @keywords internal
NULL

#' @rdname samc-class-access
#' @export
setMethod("$", signature(x = "samc"), function(x, name) {
  if(name == "override"){
    return(x@override)
  } else if (name == "q_matrix"){
    return(x@data@q)
  } else if (name == "p_matrix") {
    p <- rbind(x@data@q, rep(0, ncol(x@data@q)))
    p <- cbind(p, c(rowSums(x@data@r), 1))
    return(p)
  } else if (name == "r_matrix"){
    return(x@data@r)
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(NULL)
})


#' Modify samc-class components
#'
#' Allows users to modifu a subset of the samc-class components
#'
#' @param x samc-class object
#' @param name Component of the samc-class to modify
#' @param value Value to assign to samc-class component
#'
#' @return Updated samc-class object
#'
#' @name samc-class-modify
#' @keywords internal
NULL

#' @rdname samc-class-modify
#' @export
setMethod("$<-", signature(x = "samc"), function(x, name, value) {
  if (name == "override") {
    x@override <- value
  } else if (name == "q_matrix"){
    warning("Cannot modify the Q matrix this way.", call. = FALSE)
  } else if (name == "r_matrix"){
    warning("Cannot modify the R matrix this way.", call. = FALSE)
  } else if (name == "p_matrix"){
    warning("Cannot modify the P matrix this way.", call. = FALSE)
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(x)
})
