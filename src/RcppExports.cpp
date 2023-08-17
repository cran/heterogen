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

// bg_transpose
NumericVector bg_transpose(NumericMatrix mat, const bool byrow);
RcppExport SEXP _heterogen_bg_transpose(SEXP matSEXP, SEXP byrowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const bool >::type byrow(byrowSEXP);
    rcpp_result_gen = Rcpp::wrap(bg_transpose(mat, byrow));
    return rcpp_result_gen;
END_RCPP
}
// distance_weighted_gauss
arma::vec distance_weighted_gauss(arma::mat coord_xy, arma::vec point_xy, double tau);
RcppExport SEXP _heterogen_distance_weighted_gauss(SEXP coord_xySEXP, SEXP point_xySEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coord_xy(coord_xySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type point_xy(point_xySEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(distance_weighted_gauss(coord_xy, point_xy, tau));
    return rcpp_result_gen;
END_RCPP
}
// float_relative
NumericVector float_relative(NumericVector xx);
RcppExport SEXP _heterogen_float_relative(SEXP xxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    rcpp_result_gen = Rcpp::wrap(float_relative(xx));
    return rcpp_result_gen;
END_RCPP
}
// float_round
NumericVector float_round(const NumericVector& float_n, int digits);
RcppExport SEXP _heterogen_float_round(SEXP float_nSEXP, SEXP digitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type float_n(float_nSEXP);
    Rcpp::traits::input_parameter< int >::type digits(digitsSEXP);
    rcpp_result_gen = Rcpp::wrap(float_round(float_n, digits));
    return rcpp_result_gen;
END_RCPP
}
// matrixmult
SEXP matrixmult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _heterogen_matrixmult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixmult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matrixvec_plus
NumericMatrix matrixvec_plus(NumericMatrix& X, NumericVector& y);
RcppExport SEXP _heterogen_matrixvec_plus(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrixvec_plus(X, y));
    return rcpp_result_gen;
END_RCPP
}
// matrixvec_subs
NumericMatrix matrixvec_subs(NumericMatrix& X, NumericVector& y);
RcppExport SEXP _heterogen_matrixvec_subs(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrixvec_subs(X, y));
    return rcpp_result_gen;
END_RCPP
}
// matrixcec_square
NumericMatrix matrixcec_square(NumericMatrix& X, NumericVector& y);
RcppExport SEXP _heterogen_matrixcec_square(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matrixcec_square(X, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_heterogen_bg_transpose", (DL_FUNC) &_heterogen_bg_transpose, 2},
    {"_heterogen_distance_weighted_gauss", (DL_FUNC) &_heterogen_distance_weighted_gauss, 3},
    {"_heterogen_float_relative", (DL_FUNC) &_heterogen_float_relative, 1},
    {"_heterogen_float_round", (DL_FUNC) &_heterogen_float_round, 2},
    {"_heterogen_matrixmult", (DL_FUNC) &_heterogen_matrixmult, 2},
    {"_heterogen_matrixvec_plus", (DL_FUNC) &_heterogen_matrixvec_plus, 2},
    {"_heterogen_matrixvec_subs", (DL_FUNC) &_heterogen_matrixvec_subs, 2},
    {"_heterogen_matrixcec_square", (DL_FUNC) &_heterogen_matrixcec_square, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_heterogen(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
