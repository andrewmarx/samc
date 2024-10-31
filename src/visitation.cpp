// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <Rcpp/Benchmark/Timer.h>

#include "solver-cache.h"
#include "constants.h"

// [[Rcpp::export(".sum_qpow_row")]]
Rcpp::List sum_qpow_row(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                        const Eigen::Map<Eigen::VectorXd> &vec,
                        const Rcpp::NumericVector &t)
{
  int n = t.size();

  Rcpp::List res(n - 1);

  Eigen::RowVectorXd vecq = vec;
  Eigen::RowVectorXd time_res = vec;

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for (int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      vecq = vecq * M;
      time_res += vecq;
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}

// [[Rcpp::export(".sum_qpow_col")]]
Rcpp::List sum_qpow_col(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                        const Eigen::Map<Eigen::VectorXd> &vec,
                        const Rcpp::NumericVector &t)
{
  int n = t.size();

  Rcpp::List res(n - 1);

  Eigen::VectorXd qc = vec;
  Eigen::VectorXd time_res = vec;

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for (int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      qc = M * qc;
      time_res += qc;
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}

// [[Rcpp::export(".f_row")]]
Rcpp::NumericVector f_row(const Eigen::SparseMatrix<double> &M,
                          const Eigen::VectorXd &vec,
                          Rcpp::XPtr<SolverCache> &SC)
{
  SC->buildSolver(M.transpose(), "mt");

  Eigen::VectorXd res = SC->solver().solve(vec);
  if(SC->solver().info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f_row");
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_row_iter")]]
Rcpp::NumericVector f_row_iter(const Eigen::SparseMatrix<double> &M,
                               const Eigen::VectorXd &vec)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M.transpose());
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in f_row_iter");
  }

  Eigen::VectorXd res = solver.solve(vec);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f_row_iter");
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col")]]
Rcpp::NumericVector f_col(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                          const Eigen::VectorXd &vec,
                          Rcpp::XPtr<SolverCache> &SC)
{
  int sz = M.rows();

  SC->buildSolver(M, "m");

  Eigen::VectorXd res = SC->solver().solve(vec);
  if(SC->solver().info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f_col");
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".f_col_iter")]]
Rcpp::NumericVector f_col_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                               const Eigen::VectorXd &vec)
{
  int sz = M.rows();

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in f_col_iter");
  }

  Eigen::VectorXd res = solver.solve(vec);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solver failed in f_col_iter");
  }

  return Rcpp::wrap(res);
}
