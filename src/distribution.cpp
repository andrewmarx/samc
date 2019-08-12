// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export(".qpow_row")]]
Rcpp::NumericVector qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int row, const int steps)
{
  Eigen::RowVectorXd res = M.row(row-1);

  for(int i = 1; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();
    res = res * M;
  }

  return Rcpp::wrap(res);
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
