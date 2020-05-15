// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::export(".sum_qpow_row")]]
Rcpp::List sum_qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int row, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd vq = Eigen::RowVectorXd::Zero(M.rows());
  vq(row - 1) = 1;

  Eigen::RowVectorXd time_res = vq;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      vq = vq * M;
      time_res = time_res + vq;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".sum_qpowrv")]]
Rcpp::List sum_qpowrv(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &rv, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd qrv = rv;
  Eigen::VectorXd time_res = qrv;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      qrv = M * qrv;
      time_res = time_res + qrv;
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

// [[Rcpp::export(".psif")]]
Rcpp::NumericVector psif(Eigen::Map<Eigen::SparseMatrix<double> > &M, Eigen::VectorXd &psi)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd res = solver.solve(psi);

  return Rcpp::wrap(res);
}
