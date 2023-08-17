if (getRversion() >= "2.15.1") { utils::globalVariables(c("bg_xy","xy_raw","envin","envin_t","bandwidth"))}

#' Core Function of GWPCA
#'
#' @description The `gwpca_core` function is a core implementation of Generalized Weighted Principal Component Analysis for each iteration.
#' 
#' @param xy  A matrix containing the coordinates of the points where environmental measurements were taken. The matrix should have two columns, representing the X and Y coordinates.
#' @param p_xy A matrix containing the coordinates of the point GWPCA will be estimated. It should have two columns for X and Y coordinates.
#' @param env A data matrix representing the environmental variables. Rows represent observations (points or grid cells), and columns represent environmental variables.
#' @param env_trans Transpose of `env` matrix.
#' @param tau The bandwidth parameter for spatial weighting in GWPCA. It determines the extent of spatial influence on the estimation of principal components.
#' 
#' @return A vector of eigenvalues from local PCA
#' 
#' 
#' @export
#' 
#' 


gwpca_core <- function(xy, p_xy, env, env_trans, tau){

  wwtest <-  distance_weighted_gauss(as.matrix(xy),as.matrix(p_xy), tau)

  m1test <- float_relative(as.vector(wwtest))
  wtest <-  float_round(as.vector(m1test),10)
  gh <- env * wtest
  y.bar <- colSums(as.matrix(gh))
  z_tst <- env_trans - y.bar
  co_test <- t(z_tst)
  no_test <- co_test * wtest
  cvs_test <- matrixmult(z_tst, no_test)
  
  if(anyNA(cvs_test)){
    z <- rep(0, ncol(env))
  } else{
    z <- princomp(covmat=cvs_test)$sdev
  }
  
  return(z)
}


#  #' SpatHetero_int
#  #' @name SpatHetero_int-class
#  #' @rdname SpatHetero_int-class
#  #' @slot matrix A Matrix of EigenValues
#  #' @slot rasters A SpatRaster for each component
#  setClass("SpatHetero_int",
#           slots = c(matrix="matrix",
#                    rasters="SpatRaster"))



#' SpatHetero
#' @name SpatHetero-class
#' @rdname SpatHetero-class
#' @slot hetero A Heterogeneity Layer
#' @slot matrix SpatHetero_in data
#' @slot rasters A SpatRaster for Each Component
#' @importClassesFrom terra SpatRaster
#' @export
setClass("SpatHetero",representation(hetero = "SpatRaster", matrix = "matrix", rasters = "SpatRaster"))

 #' Plot Heterogeneity Layer
 #' 
 #' Plot 
 #' 
 #' @param x SpatHetero Class
 #' @param comp integer. Plot specific component of the heterogeneity.
 #' @param ... Plot parameters forwarded.
 #' @return No return value, called for side effects.
 #' @export
 setMethod("plot",
           signature(x = "SpatHetero"),
          function(x, comp = NULL, ...) {
            if(is.null(comp)){
              terra::plot(x@hetero)
            }
            if(!is.null(comp)){
              terra::plot(x@rasters[[comp]])
            }
          })


# #' Plot heterogeneity Layer
# #' 
# #' Plot raster
# #' 
# #' @param x SpatHetero-Class
# #' @param comp Component fo raw heterogeneity
# #' @param ... Complement parameter of plot function
# #' 
# #' @export
# plot.SpatHetero <- function(x, comp = NULL, ...) {
#   data <- x@hetero
#  if(is.null(comp)){
#    terra::plot(x@hetero)
#  }
#  if(!is.null(comp)){
#    terra::plot(x@rasters[[comp]])
#  }
# }




