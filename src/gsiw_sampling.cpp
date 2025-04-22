// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include "eta_progress_bar.hpp"

using namespace Rcpp;
using namespace arma;

template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}


// [[Rcpp::export]]
List gSIW_sampling(const arma::mat& Gamma_init, const arma::vec& eigval_H_0, const arma::mat& eigvec_H_0, 
                   const arma::vec& r_vec, const int k, const int iter, const int burnin, const int thin,
                   const bool progress = true) {
  
  
  int p = eigval_H_0.n_elem;
  arma::mat mod_Gamma_t = Gamma_init;
  List Lambda_list;
  List Gamma_list;
  
  ETAProgressBar epb;
  Progress pg(iter, progress, epb);
  for (int t = 0; t < iter; ++t) {
    if ( pg.increment() ) {
      
      // Update lambda values
      arma::vec c = arma::sum(mod_Gamma_t.each_col() % eigval_H_0 % mod_Gamma_t, 0).t() / 2;
      arma::vec Lambda_new(p);
      
      for (int i = 0; i < p; ++i) {
        Lambda_new(i) = 1 / R::rgamma(r_vec(i) - 1, 1 / c(i));
      }
      Lambda_new = sort(Lambda_new, "descend");
      
      // Update Gamma values
      arma::uvec var_idx = arma::regspace<arma::uvec>(0, p-1);
      while (!var_idx.is_empty()) {
        arma::uvec two_idx;
        if (var_idx.n_elem == 1) {
          arma::uvec remaining_idx = arma::regspace<arma::uvec>(0, p - 1);
          remaining_idx.shed_rows(var_idx); 
          two_idx = join_vert(var_idx, remaining_idx(arma::randi<arma::uvec>(1, arma::distr_param(0, remaining_idx.n_elem - 1))));
        } else {
          two_idx = arma::randperm(var_idx.n_elem, 2);
        }
        two_idx = sort(two_idx);
        
        arma::mat sub_mod_Gamma = mod_Gamma_t.rows(two_idx);
        arma::mat sub_mat = sub_mod_Gamma.each_row() / Lambda_new.t() * sub_mod_Gamma.t();
        
        arma::vec eigval_sub_mat;
        arma::mat eigvec_sub_mat;
        arma::eig_sym(eigval_sub_mat, eigvec_sub_mat, sub_mat);
        
        
        double diff_s = eigval_sub_mat(1) - eigval_sub_mat(0);
        double diff_h = eigval_H_0(two_idx(1)) - eigval_H_0(two_idx(0));
        
        double c0 = -0.5 * std::abs(diff_s * diff_h);
        
        
        double omega = std::atan(eigvec_sub_mat(1, 1)/ eigvec_sub_mat(0, 1)); //Eigenvalue 오름차순으로 되어 있어서
        
        
        double alpha0 = R::rbeta(0.5, 0.5);
        while (R::runif(0, 1) > std::exp(c0 * alpha0)) {
          alpha0 = R::rbeta(0.5, 0.5);
        }
        
        
        
        double theta = std::acos(2 * alpha0 - 1) / 2 - omega;
        if (theta < -M_PI / 2) {
          theta += M_PI;
        } else if (theta > M_PI / 2) {
          theta -= M_PI;
        }
        
        arma::vec random_signs = 2 * arma::randi<arma::vec>(2, distr_param(0, 1)) - 1;
        arma::mat rot = {{std::cos(theta), -std::sin(theta)}, 
        {std::sin(theta), std::cos(theta)}};
        mod_Gamma_t.rows(two_idx) = arma::diagmat(random_signs) * rot * mod_Gamma_t.rows(two_idx);
        
        if(var_idx.n_elem == 1){
          break;
        }
        var_idx.shed_rows(two_idx);
      }
      
      
      if (t >= burnin && t % thin == 0) {
        
        Lambda_list.push_back(arma2vec(Lambda_new));
        arma::mat Gamma_new = eigvec_H_0 * mod_Gamma_t;
        Gamma_list.push_back(Gamma_new);
      }
      
      Rcpp::checkUserInterrupt();
      
      if (Progress::check_abort()){
        return -1;
      }
        
    }
  }
   
  return List::create(Named("Lambda") = Lambda_list, 
                      Named("Gamma") = Gamma_list);
}


