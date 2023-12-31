# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' bg_transpose
#'
#' Transpose of a matrix based on row or column index.
#' 
#' @param mat A Matrix.
#' @param byrow `FALSE` computes based on row index. `TRUE` computes based on column index.
#' @return A matrix transposed.
bg_transpose <- function(mat, byrow = FALSE) {
    .Call(`_heterogen_bg_transpose`, mat, byrow)
}

#' distance_weighted_gauss
#'
#' Weighted Distance based on Gaussian function
#' 
#' @param coord_xy A Matrix with lon/lat coordinates.
#' @param point_xy lon/lat coordinate.
#' @param tau bandwidth.
#' @return A vector.
distance_weighted_gauss <- function(coord_xy, point_xy, tau) {
    .Call(`_heterogen_distance_weighted_gauss`, coord_xy, point_xy, tau)
}

#' float_relative
#' 
#' Relative sum formula
#' 
#' @param xx A Matrix with lon/lat coordinates.
#' @return A vector.
float_relative <- function(xx) {
    .Call(`_heterogen_float_relative`, xx)
}

#' float_round
#' 
#' Rouding of Numbers
#' 
#' @param float_n A numeric vector.
#' @param digits integer indicating the number of decimal places .
#' @return A vector.
float_round <- function(float_n, digits = 0L) {
    .Call(`_heterogen_float_round`, float_n, digits)
}

#' matrixmult
#'
#' Matrix Multiplication
#' 
#' @param A A Matrix.
#' @param B A Matrix
#' @return A matrix.
matrixmult <- function(A, B) {
    .Call(`_heterogen_matrixmult`, A, B)
}

#' matrixvec_plus
#'
#' Matrix Multiplication
#' 
#' @param X A Matrix.
#' @param y A Vector
#' @return A matrix.
matrixvec_plus <- function(X, y) {
    .Call(`_heterogen_matrixvec_plus`, X, y)
}

#' matrixvec_subs
#'
#' Matrix Substraction 
#' 
#' @param X A Matrix.
#' @param y A Vector
#' @return A matrix.
matrixvec_subs <- function(X, y) {
    .Call(`_heterogen_matrixvec_subs`, X, y)
}

#' matrixcec_square
#'
#' Matrix Square 
#' 
#' @param X A Matrix.
#' @param y A Vector
#' @return A matrix.
matrixcec_square <- function(X, y) {
    .Call(`_heterogen_matrixcec_square`, X, y)
}

