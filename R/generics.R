# Copyright (c) 2024 Andrew Marx. All rights reserved.
# Licensed under AGPLv3.0. See LICENSE file in the project root for details.

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
  } else if (name == "names"){
    return(x@names)
  } else if (name == "q_matrix"){
    q_mat = x@data@f

    Matrix::diag(q_mat) <- Matrix::diag(q_mat) - 1
    q_mat@x <- -q_mat@x

    rownames(q_mat) = x$names
    colnames(q_mat) = x$names
    return(q_mat)
  } else if (name == "p_matrix") {
    p <- cbind(x$q_matrix, x@data@t_abs)
    p <- rbind(p, matrix(0, nrow = 1, ncol = ncol(p)))
    Matrix::diag(p)[ncol(p)] <- 1

    if (!is.null(x$names)) {
      colnames(p) <- c(x$names, "total")
      rownames(p) <- colnames(p)
    }

    return(p)
  } else if (name == "threads"){
    return(x@threads)
  } else if (name == "solver"){
    return(x@solver)
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(NULL)
})


#' Modify samc-class components
#'
#' Allows users to modify a subset of the samc-class components
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
  } else if (name == "abs_states"){
    if (is.logical(value) && length(value) == 1 && is.na(value)) {
      x@data@c_abs <- matrix(ncol = 0, nrow = 0)
    } else {
      x@data@c_abs <- .process_abs_states(x, value)
    }
  } else if (name == "names"){
    warning("Cannot modify the transient state names.", call. = FALSE)
  } else if (name == "q_matrix"){
    warning("Cannot modify the Q matrix this way.", call. = FALSE)
  } else if (name == "p_matrix"){
    warning("Cannot modify the P matrix this way.", call. = FALSE)
  } else if (name == "threads"){
    if (is.numeric(value) && length(value) == 1 && value%%1 == 0 && value > 0) {
      x@threads <- value
      message("Important: When setting the number of threads, make sure to use a reasonable number for the machine the code will run on.\n",
              "Using the maximum number of threads supported by your hardware can make other programs non-responsive. Only do this if nothing else needs to run during the analysis.\n",
              "Specifying more threads than a machine supports in hardware will likely lead to lost performance.\n")
    } else {
      warning("Input must be a single positive integer.", call. = FALSE)
    }
  } else if (name == "solver"){
    if (x %in% c("direct", "iter")) {
      x@solver = value
    } else {
      warning("Valid inputs: \"direct\" \"iter\"", call. = FALSE)
    }
  } else {
    warning("Invalid object specified.", call. = FALSE)
  }
  return(x)
})
