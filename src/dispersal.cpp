// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

#include <Rcpp/Benchmark/Timer.h>

#include <chrono>
#include <string>
#include <iomanip>

#include "solver-cache.h"


// [[Rcpp::export(".sum_qn_q")]]
Rcpp::List sum_qn_q(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                    const Eigen::Map<Eigen::SparseMatrix<double> > &M2,
                    const Eigen::VectorXd &q,
                    Rcpp::NumericVector t)
{
  int n = t.size();

  Eigen::VectorXd q2 = q;

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M2);

  Rcpp::List res = Rcpp::List::create();

  for(int i = 1; i < n; i++) {
    for(int j = t[i - 1]; j < t[i]; j++) {
      if(j % 1000 == 0) Rcpp::checkUserInterrupt();
        q2 = M * q2;
      }

    Eigen::VectorXd time_res = solver.solve(q - q2);

    res.push_back(time_res, std::to_string((int)t[i]));
  }

  return res;
}

// [[Rcpp::export(".sum_qn_q_iter")]]
Rcpp::List sum_qn_q_iter(const Eigen::Map<Eigen::SparseMatrix<double> > &M,
                    const Eigen::Map<Eigen::SparseMatrix<double> > &M2,
                    const Eigen::VectorXd &q,
                    Rcpp::NumericVector t)
{
  int n = t.size();

  Eigen::VectorXd q2 = q;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  solver.compute(M2);

  Rcpp::List res = Rcpp::List::create();

  for(int i = 1; i < n; i++) {
    for(int j = t[i - 1]; j < t[i]; j++) {
      if(j % 1000 == 0) Rcpp::checkUserInterrupt();
      q2 = M * q2;
    }

    Eigen::VectorXd time_res = solver.solve(q - q2);

    res.push_back(time_res, std::to_string((int)t[i]));
  }

  return res;
}


// [[Rcpp::export(".diagf_par")]]
Rcpp::NumericVector diagf_par(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int threads)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";

  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  int sz = M.rows();

  Eigen::VectorXd dg(sz);

  solver.compute(M);

  Rcpp::Rcout << " Complete.\n";

  Rcpp::Rcout << "Calculating matrix inverse diagonal...\n";

  // Parallel run
  std::vector<Eigen::VectorXd> xs(threads);
  for (auto &x : xs)
    x = Eigen::VectorXd::Zero(sz);

  RcppThread::ProgressCounter progress(sz, 10);

  auto fun = [&] (unsigned int i){
    RcppThread::checkUserInterrupt();

    Eigen::VectorXd& ident = xs[(i * threads) / sz];

    ident(i) = 1;
    Eigen::VectorXd col = solver.solve(ident);
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
Rcpp::NumericVector diagf_par_iter(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int threads)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";

  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;

  int sz = M.rows();

  Eigen::VectorXd dg(sz);

  solver.compute(M);

  Rcpp::Rcout << " Complete.\n";

  Rcpp::Rcout << "Calculating matrix inverse diagonal...\n";

  // Parallel run
  std::vector<Eigen::VectorXd> xs(threads);
  for (auto &x : xs)
    x = Eigen::VectorXd::Zero(sz);

  RcppThread::ProgressCounter progress(sz, 10);

  auto fun = [&] (unsigned int i){
    RcppThread::checkUserInterrupt();

    Eigen::VectorXd& ident = xs[(i * threads) / sz];

    ident(i) = 1;
    Eigen::VectorXd col = solver.solve(ident);
    dg(i) = col(i);
    ident(i) = 0;

    progress++;
  };

  RcppThread::parallelFor(0, sz, fun, threads, threads);

  Rcpp::Rcout << "\rComplete                                                      \n";
  Rcpp::Rcout << "Diagonal has been cached. Continuing with metric calculation...\n";

  return Rcpp::wrap(dg);
}
