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


// C++ 11 compliant progress output that should be thread safe
class progressCounter {
public:
  progressCounter(size_t sz, double f) :
    lastTime(std::chrono::steady_clock::now()),
    startTime(std::chrono::steady_clock::now()),
    size(sz),
    freq(f) {}

  void operator++(int) {
    progress++;

    auto currentTime = std::chrono::steady_clock::now();
    double deltaTime = std::chrono::duration<double>(currentTime - lastTime).count();

    std::lock_guard<std::mutex> lock(mutex);
    if(deltaTime > 5.0) {
      lastTime = currentTime;

      double t = ((double)(std::chrono::duration<double>(currentTime - startTime).count() / progress) * (double)(size - progress));

      std::string units = " seconds";

      if (t > 86400) {
        t = t / 86400;
        units = " days";
      } else if (t > 3600) {
        t = t / 3600;
        units = " hours";
      } else if (t > 60) {
        t = t / 60;
        units = " minutes";
      }

      RcppThread::Rcout << "\r" << progress << "/" << size << " : " << t << units << " remaining                              ";
    }
  }

private:
  std::mutex mutex;
  std::atomic<size_t> progress{0};
  std::size_t size;
  double freq;

  std::chrono::time_point<std::chrono::steady_clock> lastTime;
  std::chrono::time_point<std::chrono::steady_clock> startTime;
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
    for(int j = t[i - 1]; j < t[i]; j++) {
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

  progressCounter progress(sz, 10.0);

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


// [[Rcpp::export(".psid_long")]]
Rcpp::NumericVector psid_long(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::VectorXd &psi, const Eigen::VectorXd &dg)
{
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  solver.compute(M.transpose());

  Eigen::VectorXd psiF = solver.solve(psi) - psi;

  return Rcpp::wrap(psiF.cwiseQuotient(dg));
}
