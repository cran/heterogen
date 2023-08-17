#define EIGEN_WARNINGS_DISABLED
#include <RcppEigen.h>

//' matrixmult
//'
//' Matrix Multiplication
//' 
//' @param A A Matrix.
//' @param B A Matrix
//' @return A matrix.
// [[Rcpp::export]]
SEXP matrixmult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}



