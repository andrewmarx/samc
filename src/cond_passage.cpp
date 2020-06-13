// Copyright (c) 2020 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::export(".cond_t")]]
Rcpp::NumericVector cond_t(Eigen::Map<Eigen::SparseMatrix<double> > &IQ, Eigen::VectorXd &qj)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(IQ);

  Eigen::VectorXd b = solver.solve(qj);

  solver.compute(IQ * b.asDiagonal());

  Eigen::VectorXd res = solver.solve(b);

  return Rcpp::wrap(res);
}


