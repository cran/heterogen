#define BOOST_DISABLE_ASSERTS

#include <Rcpp.h>
using namespace Rcpp;

//' float_round
//' 
//' Rouding of Numbers
//' 
//' @param float_n A numeric vector.
//' @param digits integer indicating the number of decimal places .
//' @return A vector.
// [[Rcpp::export]]
NumericVector float_round(const NumericVector& float_n, int digits = 0) {
  NumericVector B = clone(float_n);
  std::size_t K = float_n.size();
  for (std::size_t k = 0; k < K; k++) {
    B[k] = ::Rf_fround(float_n[k], digits);
  }
  return B;
}
