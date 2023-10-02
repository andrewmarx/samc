// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "samc_types.h"
#include <RcppEigen.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cond_t
Rcpp::NumericVector cond_t(Eigen::Map<Eigen::SparseMatrix<double> >& IQ, Eigen::VectorXd& qj);
RcppExport SEXP _samc_cond_t(SEXP IQSEXP, SEXP qjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type IQ(IQSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type qj(qjSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_t(IQ, qj));
    return rcpp_result_gen;
END_RCPP
}
// cond_t_iter
Rcpp::NumericVector cond_t_iter(Eigen::Map<Eigen::SparseMatrix<double> >& IQ, Eigen::VectorXd& qj);
RcppExport SEXP _samc_cond_t_iter(SEXP IQSEXP, SEXP qjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type IQ(IQSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type qj(qjSEXP);
    rcpp_result_gen = Rcpp::wrap(cond_t_iter(IQ, qj));
    return rcpp_result_gen;
END_RCPP
}
// build_convolution_cache
Rcpp::XPtr<convolution_cache> build_convolution_cache(const Rcpp::NumericMatrix& kernel, const Rcpp::NumericMatrix& resistance, const Rcpp::NumericMatrix& fidelity, const Rcpp::NumericMatrix& absorption, const bool symmetric, const int threads);
RcppExport SEXP _samc_build_convolution_cache(SEXP kernelSEXP, SEXP resistanceSEXP, SEXP fidelitySEXP, SEXP absorptionSEXP, SEXP symmetricSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type resistance(resistanceSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type fidelity(fidelitySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type absorption(absorptionSEXP);
    Rcpp::traits::input_parameter< const bool >::type symmetric(symmetricSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_convolution_cache(kernel, resistance, fidelity, absorption, symmetric, threads));
    return rcpp_result_gen;
END_RCPP
}
// get_convolution_list
Rcpp::List get_convolution_list(const Rcpp::XPtr<convolution_cache>& ca);
RcppExport SEXP _samc_get_convolution_list(SEXP caSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr<convolution_cache>& >::type ca(caSEXP);
    rcpp_result_gen = Rcpp::wrap(get_convolution_list(ca));
    return rcpp_result_gen;
END_RCPP
}
// convolution_short
Rcpp::List convolution_short(std::vector<long> steps, const Rcpp::XPtr<convolution_cache>& ca, const Rcpp::NumericVector& pop_in, const int threads);
RcppExport SEXP _samc_convolution_short(SEXP stepsSEXP, SEXP caSEXP, SEXP pop_inSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<long> >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::XPtr<convolution_cache>& >::type ca(caSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type pop_in(pop_inSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(convolution_short(steps, ca, pop_in, threads));
    return rcpp_result_gen;
END_RCPP
}
// convolution_long
Rcpp::List convolution_long(const Rcpp::XPtr<convolution_cache>& ca, const Rcpp::NumericVector& pop_in, const int threads);
RcppExport SEXP _samc_convolution_long(SEXP caSEXP, SEXP pop_inSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::XPtr<convolution_cache>& >::type ca(caSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type pop_in(pop_inSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(convolution_long(ca, pop_in, threads));
    return rcpp_result_gen;
END_RCPP
}
// sum_qn_q
Rcpp::List sum_qn_q(const Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::SparseMatrix<double> >& M2, const Eigen::VectorXd& q, Rcpp::NumericVector t);
RcppExport SEXP _samc_sum_qn_q(SEXP MSEXP, SEXP M2SEXP, SEXP qSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::SparseMatrix<double> >& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_qn_q(M, M2, q, t));
    return rcpp_result_gen;
END_RCPP
}
// sum_qn_q_iter
Rcpp::List sum_qn_q_iter(const Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::SparseMatrix<double> >& M2, const Eigen::VectorXd& q, Rcpp::NumericVector t);
RcppExport SEXP _samc_sum_qn_q_iter(SEXP MSEXP, SEXP M2SEXP, SEXP qSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::SparseMatrix<double> >& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_qn_q_iter(M, M2, q, t));
    return rcpp_result_gen;
END_RCPP
}
// diagf_par
Rcpp::NumericVector diagf_par(Eigen::Map<Eigen::SparseMatrix<double> >& M, const int threads);
RcppExport SEXP _samc_diagf_par(SEXP MSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(diagf_par(M, threads));
    return rcpp_result_gen;
END_RCPP
}
// diagf_par_iter
Rcpp::NumericVector diagf_par_iter(Eigen::Map<Eigen::SparseMatrix<double> >& M, const int threads);
RcppExport SEXP _samc_diagf_par_iter(SEXP MSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(diagf_par_iter(M, threads));
    return rcpp_result_gen;
END_RCPP
}
// qpow_row
Rcpp::List qpow_row(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::VectorXd>& vec, Rcpp::NumericVector steps);
RcppExport SEXP _samc_qpow_row(SEXP MSEXP, SEXP vecSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(qpow_row(M, vec, steps));
    return rcpp_result_gen;
END_RCPP
}
// qpow_col
Rcpp::List qpow_col(Eigen::Map<Eigen::SparseMatrix< double> >& M, const Eigen::Map<Eigen::VectorXd>& vec, Rcpp::NumericVector steps);
RcppExport SEXP _samc_qpow_col(SEXP MSEXP, SEXP vecSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix< double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(qpow_col(M, vec, steps));
    return rcpp_result_gen;
END_RCPP
}
// sum_qpowrv
Rcpp::List sum_qpowrv(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::VectorXd>& rv, Rcpp::NumericVector steps);
RcppExport SEXP _samc_sum_qpowrv(SEXP MSEXP, SEXP rvSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type rv(rvSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_qpowrv(M, rv, steps));
    return rcpp_result_gen;
END_RCPP
}
// solver_cache
Rcpp::XPtr<SolverCache> solver_cache();
RcppExport SEXP _samc_solver_cache() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(solver_cache());
    return rcpp_result_gen;
END_RCPP
}
// f1
Rcpp::NumericVector f1(Eigen::Map<Eigen::SparseMatrix<double> >& M, Rcpp::XPtr<SolverCache>& SC);
RcppExport SEXP _samc_f1(SEXP MSEXP, SEXP SCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SolverCache>& >::type SC(SCSEXP);
    rcpp_result_gen = Rcpp::wrap(f1(M, SC));
    return rcpp_result_gen;
END_RCPP
}
// f1_iter
Rcpp::NumericVector f1_iter(Eigen::Map<Eigen::SparseMatrix<double> >& M);
RcppExport SEXP _samc_f1_iter(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(f1_iter(M));
    return rcpp_result_gen;
END_RCPP
}
// sum_qpow_row
Rcpp::List sum_qpow_row(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::VectorXd>& vec, Rcpp::NumericVector steps);
RcppExport SEXP _samc_sum_qpow_row(SEXP MSEXP, SEXP vecSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_qpow_row(M, vec, steps));
    return rcpp_result_gen;
END_RCPP
}
// sum_qpow_col
Rcpp::List sum_qpow_col(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::Map<Eigen::VectorXd>& vec, Rcpp::NumericVector steps);
RcppExport SEXP _samc_sum_qpow_col(SEXP MSEXP, SEXP vecSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_qpow_col(M, vec, steps));
    return rcpp_result_gen;
END_RCPP
}
// f_row
Rcpp::NumericVector f_row(const Eigen::SparseMatrix<double>& M, const Eigen::VectorXd& vec, Rcpp::XPtr<SolverCache>& SC);
RcppExport SEXP _samc_f_row(SEXP MSEXP, SEXP vecSEXP, SEXP SCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SolverCache>& >::type SC(SCSEXP);
    rcpp_result_gen = Rcpp::wrap(f_row(M, vec, SC));
    return rcpp_result_gen;
END_RCPP
}
// f_row_iter
Rcpp::NumericVector f_row_iter(Eigen::SparseMatrix<double>& M, const Eigen::VectorXd& vec);
RcppExport SEXP _samc_f_row_iter(SEXP MSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(f_row_iter(M, vec));
    return rcpp_result_gen;
END_RCPP
}
// f_col
Rcpp::NumericVector f_col(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::VectorXd& vec, Rcpp::XPtr<SolverCache>& SC);
RcppExport SEXP _samc_f_col(SEXP MSEXP, SEXP vecSEXP, SEXP SCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<SolverCache>& >::type SC(SCSEXP);
    rcpp_result_gen = Rcpp::wrap(f_col(M, vec, SC));
    return rcpp_result_gen;
END_RCPP
}
// f_col_iter
Rcpp::NumericVector f_col_iter(Eigen::Map<Eigen::SparseMatrix<double> >& M, const Eigen::VectorXd& vec);
RcppExport SEXP _samc_f_col_iter(SEXP MSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::SparseMatrix<double> >& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(f_col_iter(M, vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_samc_cond_t", (DL_FUNC) &_samc_cond_t, 2},
    {"_samc_cond_t_iter", (DL_FUNC) &_samc_cond_t_iter, 2},
    {"_samc_build_convolution_cache", (DL_FUNC) &_samc_build_convolution_cache, 6},
    {"_samc_get_convolution_list", (DL_FUNC) &_samc_get_convolution_list, 1},
    {"_samc_convolution_short", (DL_FUNC) &_samc_convolution_short, 4},
    {"_samc_convolution_long", (DL_FUNC) &_samc_convolution_long, 3},
    {"_samc_sum_qn_q", (DL_FUNC) &_samc_sum_qn_q, 4},
    {"_samc_sum_qn_q_iter", (DL_FUNC) &_samc_sum_qn_q_iter, 4},
    {"_samc_diagf_par", (DL_FUNC) &_samc_diagf_par, 2},
    {"_samc_diagf_par_iter", (DL_FUNC) &_samc_diagf_par_iter, 2},
    {"_samc_qpow_row", (DL_FUNC) &_samc_qpow_row, 3},
    {"_samc_qpow_col", (DL_FUNC) &_samc_qpow_col, 3},
    {"_samc_sum_qpowrv", (DL_FUNC) &_samc_sum_qpowrv, 3},
    {"_samc_solver_cache", (DL_FUNC) &_samc_solver_cache, 0},
    {"_samc_f1", (DL_FUNC) &_samc_f1, 2},
    {"_samc_f1_iter", (DL_FUNC) &_samc_f1_iter, 1},
    {"_samc_sum_qpow_row", (DL_FUNC) &_samc_sum_qpow_row, 3},
    {"_samc_sum_qpow_col", (DL_FUNC) &_samc_sum_qpow_col, 3},
    {"_samc_f_row", (DL_FUNC) &_samc_f_row, 3},
    {"_samc_f_row_iter", (DL_FUNC) &_samc_f_row_iter, 2},
    {"_samc_f_col", (DL_FUNC) &_samc_f_col, 3},
    {"_samc_f_col_iter", (DL_FUNC) &_samc_f_col_iter, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_samc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
