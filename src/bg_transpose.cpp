#define BOOST_DISABLE_ASSERTS
#include <Rcpp.h>
using namespace Rcpp;

//' bg_transpose
//'
//' Transpose of a matrix based on row or column index.
//' 
//' @param mat A Matrix.
//' @param byrow `FALSE` computes based on row index. `TRUE` computes based on column index.
//' @return A matrix transposed.
// [[Rcpp::export]]
NumericVector bg_transpose(NumericMatrix mat, const bool byrow=false){
  if (byrow){
    mat = transpose(mat);
  }
  return NumericVector(mat);
}
