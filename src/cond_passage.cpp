// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export(".cond_t")]]
Rcpp::List cond_t(const Eigen::Map<Eigen::SparseMatrix<double> > &IQ,
                  const Eigen::VectorXd &qj)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(IQ);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in cond_t");
  }

  Eigen::VectorXd b = solver.solve(qj);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in cond_t (1)");
  }

  Eigen::VectorXd fb = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in cond_t (2)");
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("fb") = fb);

  return res;
}

// [[Rcpp::export(".cond_t_iter")]]
Rcpp::List cond_t_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &IQ,
                       const Eigen::VectorXd &qj)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(IQ);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in cond_t_iter");
  }

  Eigen::VectorXd b = solver.solve(qj);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in cond_t_iter (1)");
  }

  Eigen::VectorXd fb = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in cond_t_iter (2)");
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("fb") = fb);

  return res;
}
