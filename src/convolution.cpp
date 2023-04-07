// Copyright (c) 2023 Andrew Marx where applicable.
// Licensed under GPLv3.0. See LICENSE file in the project root for details.
//
// Derived from https://github.com/LandSciTech/LSTDConnect/blob/8927cac418454bb7910be8f832dffb9b62b34223/src/samc.cpp
// Which was licensed under GPLv3.0 at the time of use

#include <Rcpp.h>
#include "convolution.h"

// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>

template<bool SYMMETRIC>
inline void construct_cache(
    convolution_cache& ca,
    const std::vector<kernel_point_t>& kernel,
    const Rcpp::NumericMatrix& resistance,
    const Rcpp::NumericMatrix& fidelity,
    const Rcpp::NumericMatrix& absorption,
    const int threads = 1){

  ca.kernel_size = kernel.size();
  ca.nrow = resistance.nrow();
  ca.ncol = resistance.ncol();

  ca.cell_count = ca.ncol * ca.nrow;

  ca.movement_rate.clear();
  ca.movement_rate.resize(ca.kernel_size*ca.nrow*ca.ncol, {0});

  ca.absorption.assign(absorption.begin(),absorption.end());

  std::ptrdiff_t max_offset = 0;
  std::ptrdiff_t min_offset = 0;


  for(const auto& k_point : kernel){
    const std::ptrdiff_t offset = k_point.y_off + k_point.x_off*ca.nrow;
    ca.kernel.push_back(-offset);
    max_offset = std::max(max_offset, offset);
    min_offset = std::min(min_offset, offset);
  }

  ca.left_extra_cols  = std::max((ca.nrow-1-min_offset)/ca.nrow, {0});
  ca.right_extra_cols = std::max((ca.nrow-1+max_offset)/ca.nrow, {0});

  // TODO can be parallelized, but still very fast without it
  for(std::size_t x = 0; x<ca.ncol; x++){
    for(std::size_t y = 0; y<ca.nrow; y++){
      //Rcpp::Rcout << x << ", " << y << ", " << ca.death_rate[y+x*ca.nrow] << "\n";
      //bool printing = x==1 && y==1;

      double weighted_sum = 0;

      for(const auto& k_point : kernel){

        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;

        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          weighted_sum += k_point.num/(resistance[k_y + k_x*ca.nrow]+(SYMMETRIC*resistance[y + x*ca.nrow]));
          //if(printing) Rcpp::Rcout << "a:" << x << ", " << y << ", " <<  k_x << ", " <<  k_y << ", " <<  resistance[k_y + k_x*ca.nrow] << ", " <<  k_point.num << "\n";
        }
      }
      double t_scalar = 0;
      double t_fidelity = fidelity[y+x*ca.nrow];
      const double t_absorption = absorption[y+x*ca.nrow];
      if(weighted_sum == 0){
        t_scalar = 0;
        t_fidelity = 1.0-t_absorption;
      }else{
        t_scalar = (1.0-(t_fidelity+t_absorption))/weighted_sum;
      }

      const double l_scalar     = t_scalar;
      const double l_fidelity   = t_fidelity;

      //if(printing) Rcpp::Rcout << "b:" << weighted_sum << ", " << l_scalar << ", " << l_fidelity << ", " << t_absorption << "\n";

      //Rcpp::Rcout << "(" << x << "," << y << "," << scalar << ")\n";
      //can be simd
      for(std::size_t k = 0; k < ca.kernel_size; k++){
        const auto& k_point = kernel[k];
        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;

        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernel_size] = (l_scalar * k_point.num)/(resistance[k_y + k_x*ca.nrow]+(SYMMETRIC*resistance[y + x*ca.nrow])) + (k_point.x_off == 0 && k_point.y_off == 0)*l_fidelity;

          //if(printing) Rcpp::Rcout << "c:" << k_x << ", " <<  k_y << ", " <<  ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernel_size] << "\n";
        }
        //if no, leave it 0, 0 is the default value anyway
      }

    }
  }
}

// [[Rcpp::export(.build_convolution_cache)]]
Rcpp::XPtr<convolution_cache> build_convolution_cache(
    const Rcpp::NumericMatrix& kernel,
    const Rcpp::NumericMatrix& resistance,
    const Rcpp::NumericMatrix& fidelity,
    const Rcpp::NumericMatrix& absorption,
    const bool symmetric, const int threads = 1){

  std::vector<kernel_point_t> kv{};

  const std::ptrdiff_t k_nrow = kernel.nrow();
  const std::ptrdiff_t k_ncol = kernel.ncol();

  for(std::ptrdiff_t y=0; y<k_nrow; y++){
    for(std::ptrdiff_t x=0; x<k_ncol; x++){
      if(kernel[y+x*k_nrow] || (y-k_nrow/2 == 0 && x-k_ncol/2 == 0)){
        kv.push_back({y-k_nrow/2, x-k_ncol/2, kernel[y+x*k_nrow]});
        //Rcpp::Rcout << kernel[y+x*k_nrow] << "\n";
      }
    }
  }

  if(kv.size() == 0){kv.push_back({0,0,0.0});}

  convolution_cache* ca = new convolution_cache;
  if(symmetric){
    construct_cache<true>(*ca, kv, resistance, fidelity, absorption, threads);
  }else{
    construct_cache<false>(*ca, kv, resistance, fidelity, absorption, threads);
  }
  Rcpp::XPtr<convolution_cache> xp(ca, true);

  return xp;
}


inline void convolution_one_step(
    const convolution_cache& ca,
    const double* const pop_in,
    double* const pop_out,
    double* const vis,
    double &pop,
    const int threads = 1){

  const double* const p_in = pop_in;// + (ca.nrow * ca.left_extra_cols);
  double* const p_out = pop_out;// + (ca.nrow * ca.left_extra_cols);

  //std::vector<double> temp(threads);

  auto fun = [&] (unsigned int i){
    double acc = 0;
    for(std::size_t con = 0; con < ca.kernel_size; con++){
      acc += ca.movement_rate[i*ca.kernel_size+con]*p_in[i+ca.kernel[con]];
    }
    p_out[i] = acc;
    vis[i] = vis[i] + p_in[i];

    //temp[(i * threads) / ca.kernel_size] = temp[(i * threads) / ca.kernel_size] + acc;
  };

  RcppThread::parallelFor(0, ca.absorption.size(), fun, threads, threads);

  //pop = std::reduce(temp.begin(), temp.end());
}


// [[Rcpp::export(".get_convolution_list")]]
Rcpp::List get_convolution_list(const Rcpp::XPtr<convolution_cache>& ca)
{
  const convolution_cache& c = *ca;

  Rcpp::List res = Rcpp::List::create(
    Rcpp::Named("ncol") = c.ncol,
    Rcpp::Named("nrow") = c.nrow,
    Rcpp::Named("kernel_size") = c.kernel_size,
    Rcpp::Named("left_extra_cols") = c.left_extra_cols,
    Rcpp::Named("right_extra_cols") = c.right_extra_cols,
    Rcpp::Named("movement_rate") = c.movement_rate,
    Rcpp::Named("absorption") = c.absorption,
    Rcpp::Named("kernel") = c.kernel
  );

  return res;
}


// [[Rcpp::export(.convolution_short)]]
Rcpp::List convolution_short(
    std::vector<long> steps,
    const Rcpp::XPtr<convolution_cache>& ca,
    const Rcpp::NumericVector& pop_in,
    const int threads = 1){

  double pop = 0.0;

  std::vector<double> pop_a(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);
  std::vector<double> pop_b(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);
  std::vector<double> vis(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);

  std::memcpy(&pop_a[ca->nrow*ca->left_extra_cols], &pop_in[0], ca->nrow*ca->ncol*sizeof(double));

  std::vector<Rcpp::NumericVector> pops{};
  std::vector<Rcpp::NumericVector> visits{};

  // TODO replace pointer "views" with spans from C++20 when it becomes the standard for the package
  double* const p_a = &pop_a[ca->nrow*ca->left_extra_cols];
  double* const p_b = &pop_b[ca->nrow*ca->left_extra_cols];
  double* const vis_ptr = &vis[ca->nrow*ca->left_extra_cols];

  const double* p_in = p_a;
  double* p_out = p_b;

  long last_i = 0;
  for(long i : steps){
    for(long j=0; j<i-last_i; j++){
      convolution_one_step(*ca, p_in, p_out, vis_ptr, pop, threads);
      p_in = p_out;
      if(p_out == p_a){
        p_out = p_b;
      }else{
        p_out = p_a;
      }
    }

    pops.emplace_back(int(ca->nrow) * int(ca->ncol));
    visits.emplace_back(int(ca->cell_count));

    std::memcpy(&pops.back()[0], p_in, ca->nrow*ca->ncol * sizeof(double));
    std::memcpy(&visits.back()[0], vis_ptr, ca->nrow*ca->ncol * sizeof(double));

  }

  return Rcpp::List::create(Rcpp::Named("time") = steps, Rcpp::Named("dist") = pops, Rcpp::Named("vis") = visits);
}


// [[Rcpp::export(.convolution_long)]]
Rcpp::List convolution_long(
    const Rcpp::XPtr<convolution_cache>& ca,
    const Rcpp::NumericVector& pop_in,
    const int threads = 1){

  double pop = 1.0;
  double EPSILON = 0.0000000001;

  std::vector<double> pop_a(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);
  std::vector<double> pop_b(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);
  std::vector<double> vis(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);

  std::memcpy(&pop_a[ca->nrow*ca->left_extra_cols], &pop_in[0], ca->nrow*ca->ncol*sizeof(double));

  Rcpp::NumericVector pops (int(ca->nrow) * int(ca->ncol));
  Rcpp::NumericVector visits (int(ca->nrow) * int(ca->ncol));

  // TODO replace pointer "views" with spans from C++20 when it becomes the standard for the package
  double* const p_a = &pop_a[ca->nrow*ca->left_extra_cols];
  double* const p_b = &pop_b[ca->nrow*ca->left_extra_cols];
  double* const vis_ptr = &vis[ca->nrow*ca->left_extra_cols];

  const double* p_in = p_a;
  double* p_out = p_b;

  std::size_t count = 0;
  std::size_t count_limit = 1000000;

  do {
    count++;

//    convolution_one_step(*ca, p_in, p_out, pop, threads);
    convolution_one_step(*ca, p_in, p_out, vis_ptr, pop, threads);

    p_in = p_out;
    if(p_out == p_a){
      p_out = p_b;
    }else{
      p_out = p_a;
    }

    if(count % 100 == 0) {

      std::memcpy(&pops[0], p_in, ca->nrow*ca->ncol*sizeof(double));

      pop = Rcpp::sum(pops);
    }
    Rcpp::checkUserInterrupt();
  } while((pop > EPSILON) & (count < count_limit));

  if (count == count_limit) {
    Rcpp::Rcout << "\nConvolution iteration limit reached. Results may not have fully converged.\n";
  }

  std::memcpy(&visits[0], vis_ptr, ca->nrow*ca->ncol*sizeof(double));

//  std::memcpy(&visits[0], vis.data(), ca->nrow*ca->ncol*sizeof(double));

  return Rcpp::List::create(Rcpp::Named("time") = count, Rcpp::Named("dist") = pops, Rcpp::Named("vis") = visits);
}
