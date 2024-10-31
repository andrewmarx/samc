// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

#include <Rcpp/Benchmark/Timer.h>

#include <chrono>
#include <string>
#include <iomanip>

#include "solver-cache.h"
#include "constants.h"

// [[Rcpp::export(".sum_qn_q")]]
Rcpp::List sum_qn_q(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                    const Eigen::Map<Eigen::SparseMatrix<double> > &M2,
                    const Eigen::VectorXd &q,
                    const Rcpp::NumericVector &t)
{
  int n = t.size();

  Eigen::VectorXd q2 = q;

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M2);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in sum_qn_q");
  }

  Rcpp::List res(n - 1);

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for(int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      q2 = M * q2;
    }

    Eigen::VectorXd time_res = solver.solve(q - q2);
    if(solver.info() != Eigen::Success) {
      Rcpp::stop("Solver failed in sum_qn_q");
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}

// [[Rcpp::export(".sum_qn_q_iter")]]
Rcpp::List sum_qn_q_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                         const Eigen::Map<Eigen::SparseMatrix<double> > &M2,
                         const Eigen::VectorXd &q,
                         const Rcpp::NumericVector &t)
{
  int n = t.size();

  Eigen::VectorXd q2 = q;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M2);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in sum_qn_q_iter");
  }

  Rcpp::List res(n - 1);

  for(int i = 1; i < n; i++) {
    int t_start = t[i - 1];
    int t_end = t[i];

    for(int j = t_start; j < t_end; j++) {
      if(j % INTERRUPT_CHECK_INTERVAL == 0) {
        Rcpp::checkUserInterrupt();
      }

      q2 = M * q2;
    }

    Eigen::VectorXd time_res = solver.solve(q - q2);
    if(solver.info() != Eigen::Success) {
      Rcpp::stop("Solver failed in sum_qn_q_iter");
    }

    res[i - 1] = Rcpp::wrap(time_res);
  }

  return res;
}

// [[Rcpp::export(".diagf_par")]]
Rcpp::NumericVector diagf_par(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                              const int threads)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";
  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in diagf_par");
  }

  int sz = M.rows();
  Eigen::VectorXd dg(sz);

  Rcpp::Rcout << " Complete.\n";
  Rcpp::Rcout << "Calculating matrix inverse diagonal...\n";

  // Parallel run
  std::vector<Eigen::VectorXd> xs(threads, Eigen::VectorXd::Zero(sz));
  RcppThread::ProgressCounter progress(sz, 10);

  auto fun = [&] (unsigned int i){
    RcppThread::checkUserInterrupt();

    Eigen::VectorXd& ident = xs[(i * threads) / sz];
    ident(i) = 1;

    Eigen::VectorXd col = solver.solve(ident);
    if(solver.info() != Eigen::Success) {
      Rcpp::stop("Solver failed in diagf_par");
    }

    dg(i) = col(i);
    ident(i) = 0;
    progress++;
  };

  RcppThread::parallelFor(0, sz, fun, threads, threads);

  Rcpp::Rcout << "\rComplete                                                      \n";
  Rcpp::Rcout << "Diagonal has been cached. Continuing with metric calculation...\n";

  return Rcpp::wrap(dg);
}

// [[Rcpp::export(".diagf_par_iter")]]
Rcpp::NumericVector diagf_par_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                                   const int threads)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";
  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed in diagf_par_iter");
  }

  int sz = M.rows();
  Eigen::VectorXd dg(sz);

  Rcpp::Rcout << " Complete.\n";
  Rcpp::Rcout << "Calculating matrix inverse diagonal...\n";

  // Parallel run
  std::vector<Eigen::VectorXd> xs(threads, Eigen::VectorXd::Zero(sz));
  RcppThread::ProgressCounter progress(sz, 10);

  auto fun = [&] (unsigned int i){
    RcppThread::checkUserInterrupt();

    Eigen::VectorXd& ident = xs[(i * threads) / sz];
    ident(i) = 1;

    Eigen::VectorXd col = solver.solve(ident);
    if(solver.info() != Eigen::Success) {
      Rcpp::stop("Solver failed in diagf_par_iter");
    }

    dg(i) = col(i);
    ident(i) = 0;
    progress++;
  };

  RcppThread::parallelFor(0, sz, fun, threads, threads);

  Rcpp::Rcout << "\rComplete                                                      \n";
  Rcpp::Rcout << "Diagonal has been cached. Continuing with metric calculation...\n";

  return Rcpp::wrap(dg);
}
