// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat arma_crossprod(const arma::mat A, const arma::mat B) {
  arma::mat m = A.t() * B;
  return m;
}

// [[Rcpp::export]]
arma::mat arma_tcrossprod(const arma::mat A, const arma::mat B) {
  arma::mat m = A * B.t();
  return m;
}

// [[Rcpp::export]]
arma::mat arma_matmul(const arma::mat A, const arma::mat B){
  return A * B;
}

// [[Rcpp::export]]
arma::mat arma_matsum(const arma::mat A, const arma::mat B){
  return A + B;
}


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec arma_dmvnorm(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}


// [[Rcpp::export]]
List transform_cov(const int p, const List& Lambda, const List& Gamma){
  int n_iter = Lambda.size();
  List cov_list;
  for(int i = 0; i < n_iter; i++){
    arma::mat G = as<arma::mat>(Gamma[i]);
    arma::vec ld = as<arma::vec>(Lambda[i]);
    arma::mat L = arma::diagmat(ld);
    arma::mat cov = G * L * G.t();
    cov_list.push_back(cov);
  }
  return cov_list;
}


// [[Rcpp::export]]
List splitList(const List L){
  List A;
  int niter = L.size();
  NumericMatrix temp = as<NumericMatrix>(L[0]);
  int nrow = temp.nrow();
  int ncol = temp.ncol();
  for(int j = 0; j < ncol; j++){
    NumericMatrix res = no_init_matrix(nrow,niter);
    for(int i = 0; i < niter; i++){
      temp = as<NumericMatrix>(L[i]);
      res(_, i) = temp(_,j);
    }
    A.push_back(res);
  }
  return A;
}




