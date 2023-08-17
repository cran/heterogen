#include <Rcpp.h>
using namespace Rcpp;

//' matrixvec_plus
//'
//' Matrix Multiplication
//' 
//' @param X A Matrix.
//' @param y A Vector
//' @return A matrix.
// [[Rcpp::export]]
NumericMatrix matrixvec_plus(NumericMatrix & X, NumericVector & y){
  unsigned int ncol = X.ncol();
  unsigned int nrow = X.nrow();
  int counter = 0;
  for (unsigned int j=0; j<ncol; j++) {
    for (unsigned int i=0; i<nrow; i++)  {
      X[counter++] *= y[i];
    }
  }
  return X;
}

//' matrixvec_subs
//'
//' Matrix Substraction 
//' 
//' @param X A Matrix.
//' @param y A Vector
//' @return A matrix.
// [[Rcpp::export]]
NumericMatrix matrixvec_subs(NumericMatrix & X, NumericVector & y){
  unsigned int ncol = X.ncol();
  unsigned int nrow = X.nrow();
  int counter = 0;
  for (unsigned int j=0; j<ncol; j++) {
    for (unsigned int i=0; i<nrow; i++)  {
      X[counter++] -= y[i];
    }
  }
  return X;
}

//' matrixcec_square
//'
//' Matrix Square 
//' 
//' @param X A Matrix.
//' @param y A Vector
//' @return A matrix.
// [[Rcpp::export]]
NumericMatrix matrixcec_square(NumericMatrix & X, NumericVector & y){
  unsigned int ncol = X.ncol();
  unsigned int nrow = X.nrow();
  int counter = 0;
  for (unsigned int j=0; j<ncol; j++) {
    for (unsigned int i=0; i<nrow; i++)  {
      X[counter++] -= pow(y[i], 2.0);
    }
  }
  return X;
}
