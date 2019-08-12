// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::export(".f_row")]]
Rcpp::NumericVector f_row(Eigen::SparseMatrix<double> &M, const int row)
{
  int sz = M.rows();
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd row_vec = Eigen::VectorXd::Zero(sz);
  row_vec(row-1) = 1;

  Eigen::VectorXd res = solver.solve(row_vec);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col")]]
Rcpp::NumericVector f_col(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int col)
{
  int sz = M.rows();

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Eigen::VectorXd col_vec = Eigen::VectorXd::Zero(sz);
  col_vec(col-1) = 1;

  Eigen::VectorXd res = solver.solve(col_vec);

  return Rcpp::wrap(res);
}
