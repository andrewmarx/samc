// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Rcpp/Benchmark/Timer.h>


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
    for (int j = t[i - 1]; j < t[i]; j++) {
      if(j % 1000 == 0) Rcpp::checkUserInterrupt();
        q2 = M * q2;
      }

    Eigen::VectorXd time_res = solver.solve(q - q2);

    res.push_back(time_res, std::to_string((int)t[i]));
  }

  return res;
}


// [[Rcpp::export(".psid_long")]]
Rcpp::NumericVector psid_long(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::VectorXd &psi)
{
  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  int sz = M.rows();

  Eigen::VectorXd ident = Eigen::VectorXd::Zero(sz);
  Eigen::VectorXd dg(sz);
  Eigen::VectorXd col(sz);

  Rcpp::Timer timer;
  double t = 0.0;
  Rcpp::String units("unit");

  int ut = std::max(sz/10, 1000);
  int utm = (int)std::max(10.0, std::pow(10, (7.0 - std::log10((double)sz))));

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M);

  Rcpp::Rcout << " Complete.\n";

  Rcpp::Rcout << "Calculating matrix inverse diagonal...";

  for(int i = 0; i < sz; i++) {
    if(i % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }

    if(i % ut == utm) {
      t = ((double)((timer.now() - timer.origin()) / i) * (double)(sz - i)) / 1000000000.0;

      if (t > 86400) {
        t = t / 86400;
        units = " days";
      } else if (t > 3600) {
        t = t / 3600;
        units = " hours";
      } else if (t > 60) {
        t = t / 60;
        units = " minutes";
      } else {
        units = " seconds";
      }

      Rcpp::Rcout << "\rCalculating matrix inverse diagonal... " << t << units.get_cstring() << " remaining         ";
    }

    ident(i) = 1;
    col = solver.solve(ident);
    dg(i) = col(i);
    ident(i) = 0;
  }

  Rcpp::Rcout << "\rCalculating matrix inverse diagonal... Complete                                           \n";

  Rcpp::Rcout << "Performing final calculations. This may take a few minutes...";

  solver.compute(M.transpose());

  Eigen::VectorXd psiF = solver.solve(psi) - psi;

  Rcpp::Rcout << " Complete.\n";

  return Rcpp::wrap(psiF.cwiseQuotient(dg));
}
