// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"


// [[Rcpp::export(".sum_qpowrv")]]
Rcpp::List sum_qpowrv(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &rv, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd qrv = rv;
  Eigen::VectorXd time_res = qrv;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      qrv = M * qrv;
      time_res = time_res + qrv;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}
