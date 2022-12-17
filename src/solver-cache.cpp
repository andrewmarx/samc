// Copyright (c) 2022 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "solver-cache.h"

// [[Rcpp::export(".solver_cache")]]
Rcpp::XPtr<SolverCache> solver_cache()
{
  SolverCache* sc = new SolverCache();
  Rcpp::XPtr<SolverCache> ptr(sc);
  return ptr;
}

void SolverCache::buildSolver(Eigen::SparseMatrix<double> &M, const std::string& fun)
{
  if (name != fun)
  {
    m_solver.compute(M);
    name = fun;
  }
}
