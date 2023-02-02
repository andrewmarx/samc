// Copyright (c) 2022 Andrew Marx. All rights reserved.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.

#ifndef SOLVERCACHE_H
#define SOLVERCACHE_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>

using namespace Eigen;

class SolverCache {
  typedef SparseLU<SparseMatrix<double> > Solver;

  Solver m_solver;
  std::string name;

public:
  SolverCache() { name = ""; }
  void buildSolver(const Eigen::SparseMatrix<double> &M, const std::string& fun);
  Solver& solver() { return m_solver; }
};

#endif
