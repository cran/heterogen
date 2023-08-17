#' Heterogeneity (rasters) 
#'
#' @description The `heterog` function is designed to calculate environmental heterogeneity metric from a raster stack dataset. 
#' This function aids in assessing the spatial variation and diversity of environmental variables within the raster data,
#'  providing valuable insights into the heterogeneity of the study area.
#'  
#' @param datastack `SpatRaster` class. The input raster stack representing environmental variables. 
#' Each layer in the stack corresponds to a different environmental variable, and the function calculates heterogeneity based on the variability across these layers.
#' @param bandwidth The bandwidth for the spatial weighting function.
#' @param tolerance The tolerance for spatial weight computation.
#' @param nprocess (Optional) The number of iterations for calculating the principal components. Default is set to 1000.
#' @param ncores (Optional) The number of cores to be used for parallel computation. Only applicable if `parallel` is set to `TRUE`. Default is 4.
#' @param normalized (Optional) A logical value indicating whether the input data should be normalized before performing GWPCA. Default is `FALSE`, meaning the data will not be normalized.
#' Take in account that core function performs correlation analysis in order to normalize the input variables.
#' @param method The method used for GWPCA computation. It can take one of the following values. `local` Performs GWPCA locally and will save each iteration on .rds files. Recommended for large-scale data sets.
#' `inter` Uses RAM memory to . Default is `inter`.
#' @param dirds (Optional) The directory where the results will be saved in RDS format. Default is `rds`.
#' @param parallel (Optional) A logical value indicating whether to run the computation in parallel. If `TRUE`, multiple cores will be used for processing. Default is `FALSE`.
#' 
#' 
#' @useDynLib heterogen, .registration=TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom methods new
#' 
#' @return A SpatHetero object
#' 
#' \itemize{
#' \item{hetero A heterogeneity layer}
#' \item{matrix A Matrix of eigenvalues}
#' \item{rasters A complete set of heterogeneity layers for each component}
#' }
#' 
#' @examples
#' \donttest{
#' # Case 01: South
#' path <- system.file("extdata","south", package="heterogen")
#' south_rast <- terra::rast(list.files(path, full.names = TRUE, 
#' pattern = '.tif'))
#' 
#' south_het <- heterog(south_rast, parallel = TRUE, 
#' bandwidth = 0.1, tolerance = 10)
#' plot(south_het)
#' }
#' 
#' \donttest{
#' # Case 02: North
#' path <- system.file("extdata","north", package="heterogen")
#' north_rast <- terra::rast(list.files(path, full.names = TRUE, 
#' pattern = '.tif'))
#' 
#' north_het <- heterog(north_rast, parallel = TRUE,
#' bandwidth = 0.1, tolerance = 10)
#' plot(north_het)
#' }
#' 
#' 
#' 
#' @export
#' 

heterog <- function(datastack, bandwidth = 0.3, tolerance = 5, nprocess = 1000,  parallel = FALSE, 
                     ncores = 2 , normalized = FALSE, method = 'iter', dirds='rds'){
   
  if(!class(datastack)[1] == 'SpatRaster') 
    stop('datastack must be a SpatRaster')
  
  if(terra::res(datastack)[1] < 0) 
    stop('bandwidth must be a positive value and greater than spatial resolution')
  
  if(tolerance < 0 | tolerance > 100) 
    stop("Tolerance value must be within the range 0 to 100")
  
  ncaval <- parallel::detectCores()
  if(ncores > ncaval) 
    stop("The number of cores requested (ncores) exceeds the available number of cores on your computer")

  METHODS <- c("iter", "local")
  method <- match.arg(method, METHODS)
  
  df_vals <- terra::as.data.frame(datastack, xy=TRUE)
  
   
    zscores <-  gwpca_df_mc(datadf = as.matrix(df_vals), bandwidth, tolerance, nprocess, parallel, ncores,method=method,dirds)
   
  
  colnames(zscores)  <- paste0('het_com_',1:dim(zscores)[2])

  base_h <- datastack[[1]]
  cellsi <- terra::cellFromXY(base_h, df_vals[,1:2])
  
  hetrs <- list()
  for(i in 1:length(colnames(zscores))){
    base_h[cellsi] <- zscores[,i]
    names(base_h) <- colnames(zscores)[i]
    hetrs[[i]] <- base_h
  }
  
  hetr_summary <- apply(zscores, 1, mean)
  sm_h <- base_h; names(sm_h) <- 'heterogeneity'
  sm_h[cellsi] <-  scales::rescale(hetr_summary, to = c(0, 1))
  
  out <- new("SpatHetero", hetero=sm_h, matrix=zscores, rasters=do.call(c, hetrs))
  
  return(out)
}
