// Copyright (c) 2019-2023 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <Rcpp/Benchmark/Timer.h>

#include "solver-cache.h"

// [[Rcpp::export(".sum_qpow_row")]]
Rcpp::List sum_qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int row, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd row_vec = Eigen::RowVectorXd::Zero(M.rows());
  row_vec(row - 1) = 1;

  Eigen::RowVectorXd time_res = row_vec;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      row_vec = row_vec * M;
      time_res = time_res + row_vec;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".sum_qpow_col")]]
Rcpp::List sum_qpow_col(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int &col, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd col_vec = Eigen::VectorXd::Zero(M.rows());
  col_vec(col-1) = 1;

  Eigen::VectorXd qc = col_vec;
  Eigen::VectorXd time_res = col_vec;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      qc = M * qc;
      time_res = time_res + qc;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".sum_psiqpow")]]
Rcpp::List sum_psiqpow(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &psi, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd psiq = psi;
  Eigen::RowVectorXd time_res = psi;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      psiq = psiq * M;
      time_res = time_res + psiq;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

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
