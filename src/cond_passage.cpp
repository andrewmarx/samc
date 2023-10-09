// Copyright (c) 2020 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::export(".cond_t")]]
Rcpp::List cond_t(Eigen::Map<Eigen::SparseMatrix<double> > &IQ, Eigen::VectorXd &qj)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(IQ);

  Eigen::VectorXd b = solver.solve(qj);
  Eigen::VectorXd fb = solver.solve(b);

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("fb") = fb);

  return res;
}

// [[Rcpp::export(".cond_t_iter")]]
Rcpp::List cond_t_iter(Eigen::Map<Eigen::SparseMatrix<double> > &IQ, Eigen::VectorXd &qj)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(IQ);

  Eigen::VectorXd b = solver.solve(qj);
  Eigen::VectorXd fb = solver.solve(b);

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("fb") = fb);

  return res;
}
