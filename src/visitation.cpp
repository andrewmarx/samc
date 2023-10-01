// Copyright (c) 2019-2023 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <Rcpp/Benchmark/Timer.h>

#include "solver-cache.h"

// [[Rcpp::export(".sum_qpow_row")]]
Rcpp::List sum_qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &vec, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd vecq = vec;
  Eigen::RowVectorXd time_res = vec;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      vecq = vecq * M;
      time_res = time_res + vecq;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".sum_qpow_col")]]
Rcpp::List sum_qpow_col(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &vec, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd qc = vec;
  Eigen::VectorXd time_res = vec;

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


// [[Rcpp::export(".f_row")]]
Rcpp::NumericVector f_row(const Eigen::SparseMatrix<double> &M, const Eigen::VectorXd &vec, Rcpp::XPtr<SolverCache> &SC)
{
  SC->buildSolver(M.transpose(), "mt");

  Eigen::VectorXd res = SC->solver().solve(vec);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_row_iter")]]
Rcpp::NumericVector f_row_iter(Eigen::SparseMatrix<double> &M, const Eigen::VectorXd &vec)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd res = solver.solve(vec);

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col")]]
Rcpp::NumericVector f_col(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::VectorXd &vec, Rcpp::XPtr<SolverCache> &SC)
{
  int sz = M.rows();

  //Rcpp::Timer timer;

  //timer.step("compute() start");
  SC->buildSolver(M, "m");
  //timer.step("compute() end");

  //timer.step("solve() start");
  Eigen::VectorXd res = SC->solver().solve(vec);
  //timer.step("solve() end");

  //Rcpp::NumericVector tr(timer);
  //Rcpp::Rcout << tr;

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col_iter")]]
Rcpp::NumericVector f_col_iter(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::VectorXd &vec)
{
  int sz = M.rows();

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  //Rcpp::Timer timer;

  //timer.step("compute() start");
  solver.compute(M);
  //timer.step("compute() end");

  //timer.step("solve() start");
  Eigen::VectorXd res = solver.solve(vec);
  //timer.step("solve() end");

  //Rcpp::NumericVector tr(timer);
  //Rcpp::Rcout << tr;

  return Rcpp::wrap(res);
}


