#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


using namespace Rcpp;

// [[Rcpp::export]]
List eig_decomp(const Eigen::MatrixXd A){
  
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigen decomposition failed." << std::endl;
    return -1;
  }
  
  Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
  Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
  
  return List::create(Rcpp::Named("values") = eigenvalues, Rcpp::Named("vectors") = eigenvectors);
}

// [[Rcpp::export]]
bool is_pd(Eigen::MatrixXd A) {
  Eigen::LLT<Eigen::MatrixXd> lltOfA(A);
  return (lltOfA.info() == Eigen::Success);
}