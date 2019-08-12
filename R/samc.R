# Copyright (c) 2019 Andrew Marx and contributors. All rights reserved.
# Licensed under the MIT license. See LICENSE file in the project root for details.

#' Create an absorbing Markov chain matrix
#'
#' @param resistance A RasterLayer object.
#' @param absorption A RasterLayer object.
#' @return A spatial absorbing Markov chain matrix.
#' @export

samc <- function(resistance = NA, absorption = NA, latlon = NA, tr_fun = function(x) 1/mean(x)) {

  # Check the inputs
  if (!is(resistance, "RasterLayer")) {
    stop("The resistance parameter is not a RasterLayer object")
  }

  if (!is(absorption, "RasterLayer")) {
    stop("The absorption parameter is not a RasterLayer object")
  }

  if (!(latlon %in% c(TRUE, FALSE))) {
    stop("The latlon parameter must be set to TRUE or FALSE")
  }

  # TODO check rasters have same dimensions

  #cond <- 1 / resistance
  tr <- gdistance::transition(resistance, transitionFunction = tr_fun, 8)
  if (latlon) {
    tr <- gdistance::geoCorrection(tr, type = "c")
  } else {
    tr <- gdistance::geoCorrection(tr)
  }
  tr_mat <- gdistance::transitionMatrix(tr)
  diag(tr_mat) <- 0
  abs_vec <- as.vector(absorption)
  tr_mat <- as(tr_mat, "dgTMatrix")
  tr_mat@x <- (1 - abs_vec[tr_mat@i + 1]) * tr_mat@x / Matrix::rowSums(tr_mat)[tr_mat@i + 1]

  # Combine the transition matrix with the absorbing data
  amc_df <- data.frame(i = c(tr_mat@i, (0:(length(abs_vec) - 1))[is.finite(abs_vec)], length(abs_vec)),
                       j = c(tr_mat@j, rep(length(abs_vec), sum(is.finite(abs_vec))), length(abs_vec)),
                       x = c(tr_mat@x, abs_vec[is.finite(abs_vec)], 1))

  # 'Remove' the null values by changing the index values using a lookup vector
  lookup_vec <- 0:(length(unique(tr_mat@i)))
  names(lookup_vec) <- c(sort(unique(tr_mat@i)), length(abs_vec))

  amc_df$i <- lookup_vec[as.character(amc_df$i)]
  amc_df$j <- lookup_vec[as.character(amc_df$j)]

  # Create the final sparse matrix
  amc_mat <- list(p = Matrix::sparseMatrix(i = amc_df$i,
                                           j = amc_df$j,
                                           x = amc_df$x,
                                           index1 = FALSE))

  class(amc_mat) <- "samc"

  return(amc_mat)
}


#' Calculate the fundamental matrix
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @return Fundamental matrix.
#' @export

fundamental <- function(amc){
  check_samc(amc)

  q <- amc$p[-nrow(amc$p), -nrow(amc$p)]
  q@x <- -q@x
  diag(q) <- 1
  n <- solve(q)
  return(n)
}


#' Calculate D
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @return D matrix
#' @export

calc_d <- function(amc){
  check_samc(amc)

  f <- fundamental(amc)
  gc()
  fdg <- 1 / diag(f)
  fdg_mat <- Matrix::sparseMatrix(i = 1:length(fdg),
                          j = 1:length(fdg),
                          x = fdg,
                          index1 = TRUE)
  # TODO Check if 'diag(f) <- ' resulting in extra memory allocation.
  diag(f) <- diag(f) - 1
  gc()
  d_mat <- f %*% fdg_mat

  return(d_mat)
}


#' Calculate B
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @return B matrix
#' @export

calc_b <- function(amc){
  check_samc(amc)

  f <- fundamental(amc)
  gc()
  rdg <- amc$p[-nrow(amc$p), ncol(amc$p)]
  r <- Matrix::sparseMatrix(i = 1:length(rdg),
                    j = 1:length(rdg),
                    x = rdg,
                    index1 = TRUE)

  # TODO f %*% r can be simplified to an elementwise multiplication of the matrix columns by the corresponding elements in the rdg vector. This might be helpful for memory allocations and performance.
  b <- f %*% r
  gc()
  return(b)
}


#' Calculate transient B
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @param t
#' @return B matrix
#' @export

calc_trans_b <- function(amc, t) {
  check_samc(amc)

  if (!is.numeric(t)) {
    stop("The t parameter must be a non-negative integer")
  } else if (t %% 1 != 0) {
    stop("The t parameter must be a non-negative integer")
  } else if (t < 0) {
    stop("The t parameter must be a non-negative integer")
  }

  # TODO: remove as.matrix call, which is needed to convert from a sparse to
  # dense matrix for the %^% operator, which means removing expm as a dependency
  Q <- as.matrix(amc$p[-nrow(amc$p), -nrow(amc$p)])
  Rdiag <- amc$p[-nrow(amc$p), ncol(amc$p)]
  R <- matrix(0, nrow = length(Rdiag), ncol = length(Rdiag))
  diag(R) <- Rdiag
  I <- diag(dim(Q)[2])

  #for t-1; so if t = 3, applies to t=2
  t <- t + 1
  Q_n <- solve(I - Q) %*% (I - Q %^% t)

  Bt <- Q_n %*% R
  return(Bt)
}


#' Calculate z
#'
#' Caculates z directly without the fundamental matrix, substantially reducing the memory requirements. 1,000,000 nodes should be feasible with 8GB of ram
#'
#' @param a Absorbing markov chain matrix. This should be output from the samc() function.
#' @return z
#' @export

calc_z <- function(amc) {
  check_samc(amc)

  q = amc$p[-nrow(amc$p),-nrow(amc$p)]
  q@x <- -q@x
  diag(q) <- 1
  b <- rep(1,dim(q)[2])
  z = solve(q, b)

  return(z)
}


#' Calculate mortatility
#'
#' Caculate the product of the patch vector and B
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @param pv Patch raster or vector
#' @return mort
#' @export

calc_mort <- function(amc, pv) {
  # TODO Add parameter check for pv

  check_samc(amc)

  pv <- as.vector(pv)
  pv <- pv[is.finite(pv)]
  q = amc$p[-nrow(amc$p),-nrow(amc$p)]
  rdg <- amc$p[-nrow(amc$p), ncol(amc$p)]
  rdg <- 1 / rdg
  r <- Matrix::sparseMatrix(i = 1:length(rdg),
                            j = 1:length(rdg),
                            x = rdg,
                            index1 = TRUE)

  A <- Matrix::t(r - r %*% q)

  mort <- solve(A, pv)

  return(mort)
}


#' Calculate dispersal
#'
#' Caculate the product of the patch vector and D
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @param pv Patch raster or vector
#' @return disp
#' @export

calc_disp <- function(amc, pv) {
  # TODO Add parameter check for pv

  check_samc(amc)

  pv <- as.vector(pv)
  pv <- pv[is.finite(pv)]

  dm <- calc_d(amc)

  disp <- pv %*% dm

  return(disp)
}


#' Calculate transient mortatility
#'
#' Caculate the product of the patch vector and transient B
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @param t Time
#' @param pv Patch raster or vector
#' @return mort
#' @export

calc_trans_mort <- function(amc, t, pv) {
  # TODO Add parameter check for pv

  check_samc(amc)

  if (!is.numeric(t)) {
    stop("The t parameter must be a non-negative integer")
  } else if (t %% 1 != 0) {
    stop("The t parameter must be a non-negative integer")
  } else if (t < 0) {
    stop("The t parameter must be a non-negative integer")
  }

  pv <- as.vector(pv)
  pv <- pv[is.finite(pv)]

  tb <- calc_trans_b(amc, t)
  tbp <- pv %*% tb

  return(tbp)
}


#' Check if input is of class 'samc'
#'
#' @param amc Absorbing markov chain matrix. This should be output from the samc() function.
#' @return None

check_samc <- function(amc) {
  if (!is(amc, "samc")) stop("The amc parameter must be output from the samc() function")
}
