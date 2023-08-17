#include "heterogen.h"
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;
using namespace stats;

// model implemented

// global model

// add distance weighted functions

// add Exponential function

// add Box-car

// add Bi-square

// add Tri-cube

// add Gaussian


double sp_dists_pp(double lon1, double lon2, double lat1, double lat2, double tau) {
  double as1 = (((lon2 - lon1 ) * (lon2 - lon1 )) + ((lat2 - lat1 ) * (lat2 - lat1 )));
  double as2 = std::exp(-as1 / (2 * (tau * tau)));
  return as2;
}


//' distance_weighted_gauss
//'
//' Weighted Distance based on Gaussian function
//' 
//' @param coord_xy A Matrix with lon/lat coordinates.
//' @param point_xy lon/lat coordinate.
//' @param tau bandwidth.
//' @return A vector.
// [[Rcpp::export]]
arma::vec distance_weighted_gauss(arma::mat coord_xy, arma::vec point_xy, double tau) {
  int N = coord_xy.n_rows, j;
  vec dists(N, fill::zeros);
  double uout = point_xy(0), vout = point_xy(1);
  for (j = 0; j < N; j++) {
    
    dists(j) = sp_dists_pp(coord_xy(j, 0), uout, coord_xy(j, 1), vout, tau);
    
  }
  return dists;
}

