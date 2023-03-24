// Copyright (c) 2023 Andrew Marx where applicable.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.
//
// Derived from https://github.com/LandSciTech/LSTDConnect/blob/8927cac418454bb7910be8f832dffb9b62b34223/inst/include/samc.h
// Which was licensed under GPLv3.0 at the time of use

#ifndef __CONVOLUTION_H__
#define __CONVOLUTION_H__

#include <Rcpp.h>

struct kernel_point_t{
  std::ptrdiff_t x_off;
  std::ptrdiff_t y_off;
  double num = 1;
};

struct convolution_cache{
  std::size_t ncol;
  std::size_t nrow;
  std::size_t cell_count;
  std::size_t kernel_size;
  std::size_t left_extra_cols;
  std::size_t right_extra_cols;
  std::vector<double> movement_rate;  /* size should be ncol*nrow*kernel_size */
  std::vector<double> absorption;     /* size should be ncol*nrow or empty*/
  std::vector<std::ptrdiff_t> kernel; /* a list of offsets from a given point to all of the places that it needs to look in the list of points */

};

#endif /* __CONVOLUTION_H__ */
