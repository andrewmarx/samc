// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"

// [[Rcpp::export(".solver_cache")]]
Rcpp::XPtr<SolverCache> solver_cache()
{
  SolverCache* sc = new SolverCache();
  Rcpp::XPtr<SolverCache> ptr(sc, true);
  return ptr;
}

void SolverCache::buildSolver(const Eigen::SparseMatrix<double> &M,
                              const std::string &fun)
{
  if (name != fun) {
    m_solver.compute(M);
    if (m_solver.info() != Eigen::Success) {
      Rcpp::stop("Solver failed in buildSolver");
    }

    name = fun;
  }
}
