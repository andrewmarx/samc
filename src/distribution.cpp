// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "constants.h"

// [[Rcpp::export(".qpow_row")]]
Rcpp::List qpow_row(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                    const Eigen::Map<Eigen::VectorXd> &vec,
                    const Rcpp::NumericVector &t)
{
  int n = t.size();

  Rcpp::List res(n - 1);

  Eigen::RowVectorXd time_res = vec;

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for (int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      time_res = time_res * M;
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}

// [[Rcpp::export(".qpow_col")]]
Rcpp::List qpow_col(const Eigen::Map<Eigen::SparseMatrix< double> > &M,
                    const Eigen::Map<Eigen::VectorXd> &vec,
                    const Rcpp::NumericVector &t)
{
  int n = t.size();

  Rcpp::List res(n - 1);

  Eigen::VectorXd time_res = M * vec;

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for (int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      time_res = M * time_res;
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}
