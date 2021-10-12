// Copyright (c) 2019 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

#include <Rcpp/Benchmark/Timer.h>

#include <chrono>
#include <string>

using namespace RcppParallel;


// C++ 11 compliant stopwatch that should be thread safe
class timer {
public:
  std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
  timer() : lastTime(std::chrono::high_resolution_clock::now()) {}
  inline double elapsed() {
    std::chrono::time_point<std::chrono::high_resolution_clock> thisTime=std::chrono::high_resolution_clock::now();
    double deltaTime = std::chrono::duration<double>(thisTime - lastTime).count();
    //lastTime = thisTime;
    return deltaTime;
  }
};



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

  // Limit how freq the time updates
  int ut = std::max(sz/10, 1000);

  // Used to given plenty of time for the first time update to be reasonable
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
  // Progress update related stuff
  tbb::mutex mtx;
  int total;
  int current;

  // Setup the solver
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  // Output vector
  RVector<double> output;

  // Start time
  std::chrono::time_point<std::chrono::steady_clock> startTime;

  // Initialize
  DiagWorker(Rcpp::NumericVector output)
    : output(output) {}

  // Calculate the diagonal
  void operator()(std::size_t begin, std::size_t end) {

    mtx.lock();
    RcppThread::Rcout << "Thread " << (int)begin << " " << (int)end << "\n";
    mtx.unlock();

    Eigen::VectorXd ident = Eigen::VectorXd::Zero(output.length());
    Eigen::VectorXd col(output.length());

    float total_time = float(end-begin);

    timer stopwatch;
    double t = 0.0;

    std::string units = "";
    int ut = std::max((int)end/10, 1000);

    for(int i = begin; i < end; i++) {
      //if(i % 100 == 0) {
        RcppThread::checkUserInterrupt();
      //}

      // begin == 0 ensures only one thread does the interrupt check and progress output
      if(begin == 0) {


        if(i % ut == ut - 1) {
          ut = 100;
          t = ((double)(stopwatch.elapsed()/ i) * (double)(end - i));

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

          //RcppThread::Rcout << "\rCalculating matrix inverse diagonal... " << t << units << " remaining         ";
        }
      }

      ident(i) = 1;
      col = solver.solve(ident);

      output[i] = col(i);
      ident(i) = 0;
    }
  }
};


// [[Rcpp::export(".diagf_par")]]
Rcpp::NumericVector diagf_par(Eigen::Map<Eigen::SparseMatrix<double> > &M)
{
  Rcpp::Rcout << "\nCached diagonal not found.\n";

  Rcpp::Rcout << "Performing setup. This can take several minutes...";

  int sz = M.rows();

  Rcpp::NumericVector dg(sz);

  // Create the worker
  DiagWorker diagWorker(dg);
  diagWorker.solver.compute(M);

  Rcpp::Rcout << " Complete.\n";

  diagWorker.startTime = std::chrono::steady_clock::now();

  //Rcpp::Rcout << "Calculating matrix inverse diagonal...";

  // Parallel run
  parallelFor(0, sz, diagWorker);

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
