// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::export(".sum_qpow_row")]]
Rcpp::NumericVector sum_qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int row, const int steps)
{
  Eigen::RowVectorXd vq = Eigen::RowVectorXd::Zero(M.rows());
  vq(row - 1) = 1;

  Eigen::RowVectorXd res = vq;

  for(int i = 1; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();

    vq = vq * M;
    res = res + vq;
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".sum_qpowrv")]]
Rcpp::NumericVector sum_qpowrv(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &rv, const int steps)
{
  Eigen::VectorXd qrv = rv;
  Eigen::VectorXd res = qrv;

  for(int i = 1; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();

    qrv = M * qrv;
    res = res + qrv;
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".sum_psiqpow")]]
Rcpp::NumericVector sum_psiqpow(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &psi, const int steps)
{
  Eigen::RowVectorXd psiq = psi;
  Eigen::RowVectorXd res = psi;

  for(int i = 1; i < steps; i++) {
    if(i % 1000 == 0) Rcpp::checkUserInterrupt();

    psiq = psiq * M;
    res = res + psiq;
  }

  return Rcpp::wrap(res);
}

// [[Rcpp::export(".psif")]]
Rcpp::NumericVector psif(Eigen::Map<Eigen::SparseMatrix<double> > &M, Eigen::VectorXd &psi)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd res = solver.solve(psi);

  return Rcpp::wrap(res);
}
