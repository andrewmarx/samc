// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"


// [[Rcpp::export(".f1")]]
Rcpp::NumericVector f1(Eigen::Map<Eigen::SparseMatrix<double> > &M, Rcpp::XPtr<SolverCache> &SC)
{
  Eigen::VectorXd one(M.rows());
  one.fill(1.0);

  SC->buildSolver(M, "m");

  Eigen::VectorXd res = SC->solver().solve(one);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f1_iter")]]
Rcpp::NumericVector f1_iter(Eigen::Map<Eigen::SparseMatrix<double> > &M)
{
  Eigen::VectorXd one(M.rows());
  one.fill(1.0);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M);

  Eigen::VectorXd res = solver.solve(one);

  return Rcpp::wrap(res);
}
