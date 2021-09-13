// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include <Rcpp/Benchmark/Timer.h>

using namespace RcppParallel;



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


// [[Rcpp::export(".diagf")]]
Rcpp::NumericVector diagf(Eigen::Map<Eigen::SparseMatrix<double> > &M)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";

  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  int sz = M.rows();

  Eigen::VectorXd ident = Eigen::VectorXd::Zero(sz);
  Eigen::VectorXd dg(sz);
  Eigen::VectorXd col(sz);

  Rcpp::Timer timer;
  double t = 0.0;
  Rcpp::String units("unit");

  int ut = std::max(sz/10, 1000);
  int utm = (int)std::max(10.0, std::pow(10, (7.0 - std::log10((double)sz))));

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
  Rcpp::Rcout << "Diagonal has been cached. Continuing with metric calculation...\n";

  return Rcpp::wrap(dg);
}


struct DiagWorker : public Worker
{

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  // Output vector
  RVector<double> output;

  // initialize with source and destination
  DiagWorker(Rcpp::NumericVector output)
    : output(output) {}

  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {

    Eigen::VectorXd ident = Eigen::VectorXd::Zero(output.length());
    Eigen::VectorXd col(output.length());

    for(int i = begin; i < end; i++) {
      ident(i) = 1;
      col = solver.solve(ident);

      output[i] = col(i);
      ident(i) = 0;
    }
  }
};


// [[Rcpp::export(".diagf_par")]]
Rcpp::NumericVector diagf_par(Eigen::Map<Eigen::SparseMatrix<double> > &M, const int grain)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";

  //Rcpp::Rcout << "Performing setup. This can take several minutes...";

  //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  int sz = M.rows();

  Rcpp::NumericVector dg(sz);

  /*Rcpp::Timer timer;
  double t = 0.0;
  Rcpp::String units("unit");

  int ut = std::max(sz/10, 1000);
  int utm = (int)std::max(10.0, std::pow(10, (7.0 - std::log10((double)sz))));
*/
  //solver.compute(M);

  //Rcpp::Rcout << " Complete.\n";

  Rcpp::Rcout << "Calculating matrix inverse diagonal...";

  // Create the worker
  //DiagWorker diagWorker(solver, dg);
  //DiagWorker::solver.compute(M);
  DiagWorker diagWorker(dg);
  diagWorker.solver.compute(M);

  // Parallel run
  parallelFor(0, sz, diagWorker, grain);


  Rcpp::Rcout << "\rCalculating matrix inverse diagonal... Complete                                           \n";
  Rcpp::Rcout << "Diagonal has been cached. Continuing with metric calculation...\n";

  return Rcpp::wrap(dg);
}


// [[Rcpp::export(".psid_long")]]
Rcpp::NumericVector psid_long(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::VectorXd &psi, const Eigen::VectorXd &dg)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd psiF = solver.solve(psi) - psi;

  return Rcpp::wrap(psiF.cwiseQuotient(dg));
}
