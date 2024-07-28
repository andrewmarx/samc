// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"

// [[Rcpp::export(".f1")]]
Rcpp::NumericVector f1(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                       Rcpp::XPtr<SolverCache> &SC)
{
  Eigen::VectorXd one(M.rows());
  one.fill(1.0);

  SC->buildSolver(M, "m");

  Eigen::VectorXd res = SC->solver().solve(one);
  if(SC->solver().info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f1");
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f1_iter")]]
Rcpp::NumericVector f1_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &M)
{
  Eigen::VectorXd one(M.rows());
  one.fill(1.0);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in f1_iter");
  }

  Eigen::VectorXd res = solver.solve(one);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f1_iter");
  }

  return Rcpp::wrap(res);
}
