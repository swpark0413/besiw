// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_crossprod
arma::mat arma_crossprod(const arma::mat A, const arma::mat B);
RcppExport SEXP _besiw_arma_crossprod(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_crossprod(A, B));
    return rcpp_result_gen;
END_RCPP
}
// arma_tcrossprod
arma::mat arma_tcrossprod(const arma::mat A, const arma::mat B);
RcppExport SEXP _besiw_arma_tcrossprod(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_tcrossprod(A, B));
    return rcpp_result_gen;
END_RCPP
}
// arma_matmul
arma::mat arma_matmul(const arma::mat A, const arma::mat B);
RcppExport SEXP _besiw_arma_matmul(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_matmul(A, B));
    return rcpp_result_gen;
END_RCPP
}
// arma_matsum
arma::mat arma_matsum(const arma::mat A, const arma::mat B);
RcppExport SEXP _besiw_arma_matsum(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_matsum(A, B));
    return rcpp_result_gen;
END_RCPP
}
// arma_dmvnorm
arma::vec arma_dmvnorm(arma::mat const& x, arma::rowvec const& mean, arma::mat const& sigma, bool const logd);
RcppExport SEXP _besiw_arma_dmvnorm(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec const& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_dmvnorm(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// transform_cov
List transform_cov(const int p, const List& Lambda, const List& Gamma);
RcppExport SEXP _besiw_transform_cov(SEXP pSEXP, SEXP LambdaSEXP, SEXP GammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const List& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const List& >::type Gamma(GammaSEXP);
    rcpp_result_gen = Rcpp::wrap(transform_cov(p, Lambda, Gamma));
    return rcpp_result_gen;
END_RCPP
}
// splitList
List splitList(const List L);
RcppExport SEXP _besiw_splitList(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(splitList(L));
    return rcpp_result_gen;
END_RCPP
}
// eig_decomp
List eig_decomp(const Eigen::MatrixXd A);
RcppExport SEXP _besiw_eig_decomp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eig_decomp(A));
    return rcpp_result_gen;
END_RCPP
}
// is_pd
bool is_pd(Eigen::MatrixXd A);
RcppExport SEXP _besiw_is_pd(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(is_pd(A));
    return rcpp_result_gen;
END_RCPP
}
// gSIW_sampling
List gSIW_sampling(const arma::mat& Gamma_init, const arma::vec& eigval_H_0, const arma::mat& eigvec_H_0, const arma::vec& r_vec, const int k, const int iter, const int burnin, const int thin, const bool progress);
RcppExport SEXP _besiw_gSIW_sampling(SEXP Gamma_initSEXP, SEXP eigval_H_0SEXP, SEXP eigvec_H_0SEXP, SEXP r_vecSEXP, SEXP kSEXP, SEXP iterSEXP, SEXP burninSEXP, SEXP thinSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Gamma_init(Gamma_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eigval_H_0(eigval_H_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type eigvec_H_0(eigvec_H_0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type r_vec(r_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(gSIW_sampling(Gamma_init, eigval_H_0, eigvec_H_0, r_vec, k, iter, burnin, thin, progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_besiw_arma_crossprod", (DL_FUNC) &_besiw_arma_crossprod, 2},
    {"_besiw_arma_tcrossprod", (DL_FUNC) &_besiw_arma_tcrossprod, 2},
    {"_besiw_arma_matmul", (DL_FUNC) &_besiw_arma_matmul, 2},
    {"_besiw_arma_matsum", (DL_FUNC) &_besiw_arma_matsum, 2},
    {"_besiw_arma_dmvnorm", (DL_FUNC) &_besiw_arma_dmvnorm, 4},
    {"_besiw_transform_cov", (DL_FUNC) &_besiw_transform_cov, 3},
    {"_besiw_splitList", (DL_FUNC) &_besiw_splitList, 1},
    {"_besiw_eig_decomp", (DL_FUNC) &_besiw_eig_decomp, 1},
    {"_besiw_is_pd", (DL_FUNC) &_besiw_is_pd, 1},
    {"_besiw_gSIW_sampling", (DL_FUNC) &_besiw_gSIW_sampling, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_besiw(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
