// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <Rcpp/Benchmark/Timer.h>

#include "solver-cache.h"


// [[Rcpp::export(".f_row")]]
Rcpp::NumericVector f_row(const Eigen::SparseMatrix<double> &M, const int row, Rcpp::XPtr<SolverCache> &SC)
{
  int sz = M.rows();

  SC->buildSolver(M.transpose(), "mt");

  Eigen::VectorXd row_vec = Eigen::VectorXd::Zero(sz);
  row_vec(row-1) = 1;

  Eigen::VectorXd res = SC->solver().solve(row_vec);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_row_iter")]]
Rcpp::NumericVector f_row_iter(Eigen::SparseMatrix<double> &M, const int row)
{
  int sz = M.rows();
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd row_vec = Eigen::VectorXd::Zero(sz);
  row_vec(row-1) = 1;

  Eigen::VectorXd res = solver.solve(row_vec);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col")]]
Rcpp::NumericVector f_col(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int col, Rcpp::XPtr<SolverCache> &SC)
{
  int sz = M.rows();

  //Rcpp::Timer timer;

  //timer.step("compute() start");
  SC->buildSolver(M, "m");
  //timer.step("compute() end");

  Eigen::VectorXd col_vec = Eigen::VectorXd::Zero(sz);
  col_vec(col-1) = 1;

  //timer.step("solve() start");
  Eigen::VectorXd res = SC->solver().solve(col_vec);
  //timer.step("solve() end");

  //Rcpp::NumericVector tr(timer);
  //Rcpp::Rcout << tr;

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col_iter")]]
Rcpp::NumericVector f_col_iter(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int col)
{
  int sz = M.rows();

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  //Rcpp::Timer timer;

  //timer.step("compute() start");
  solver.compute(M);
  //timer.step("compute() end");

  Eigen::VectorXd col_vec = Eigen::VectorXd::Zero(sz);
  col_vec(col-1) = 1;

  //timer.step("solve() start");
  Eigen::VectorXd res = solver.solve(col_vec);
  //timer.step("solve() end");

  //Rcpp::NumericVector tr(timer);
  //Rcpp::Rcout << tr;

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".psif")]]
Rcpp::NumericVector psif(Eigen::Map<Eigen::SparseMatrix<double> > &M, Eigen::VectorXd &psi, Rcpp::XPtr<SolverCache> &SC)
{
  SC->buildSolver(M.transpose(), "mt");

  Eigen::VectorXd res = SC->solver().solve(psi);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".psif_iter")]]
Rcpp::NumericVector psif_iter(Eigen::Map<Eigen::SparseMatrix<double> > &M, Eigen::VectorXd &psi)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd res = solver.solve(psi);

  return Rcpp::wrap(res);
}
