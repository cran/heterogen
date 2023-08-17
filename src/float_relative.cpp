#include <RcppEigen.h>
//#include "heterogen_types.h"

using namespace Rcpp;

using  Eigen::Map;
using  Eigen::VectorXd;
typedef  Map<VectorXd>  MapVecd;

// [[Rcpp::plugins(cpp11)]]

//' float_relative
//' 
//' Relative sum formula
//' 
//' @param xx A Matrix with lon/lat coordinates.
//' @return A vector.
// [[Rcpp::export]]
NumericVector float_relative(NumericVector xx) {
  const MapVecd x(as<MapVecd>(xx));
  return xx/x.sum();
}
