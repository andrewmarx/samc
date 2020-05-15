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
Rcpp::List qpow_col(Eigen::Map<Eigen::SparseMatrix< double> > &M, const int col, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd time_res = M.col(col-1);

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      time_res = M * time_res;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".psiq")]]
Rcpp::List psiq(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &psi, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd time_res = psi;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      time_res = time_res * M;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}
