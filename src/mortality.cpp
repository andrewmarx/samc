// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"
#include "constants.h"

// [[Rcpp::export(".sum_qpowrv")]]
Rcpp::List sum_qpowrv(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                      const Eigen::Map<Eigen::VectorXd> &rv,
                      const Rcpp::NumericVector &t)
{
  int n = t.size();

  Rcpp::List res(n - 1);

  Eigen::VectorXd qrv = rv;
  Eigen::VectorXd time_res = qrv;

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for (int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      qrv = M * qrv;
      time_res += qrv;
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}
