// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export(".qpow_row")]]
Rcpp::List qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int row, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd time_res = M.row(row-1);

  for(int i = 1; i < n; i++) {

    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      time_res = time_res * M;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".qpow_col")]]
Rcpp::NumericVector qpow_col(Eigen::Map<Eigen::SparseMatrix< double> > &M, const int col, const int steps)
{
  Eigen::VectorXd res = M.col(col-1);

  for(int i = 1; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();
    res = M * res;
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".psiq")]]
Rcpp::NumericVector psiq(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &psi, const int steps)
{
  Eigen::RowVectorXd res = psi;

  for(int i = 0; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();
    res = res * M;
  }

  return Rcpp::wrap(res);
}
