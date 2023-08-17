#' Perform GWPCA from data.frame with spatial structure. 
#'
#' @description gwpca_df is an R function that performs Generalized Weighted Principal Component Analysis (GWPCA) on a given dataset. 
#' This function allow to calculate the environmental heterogeneity from data.frame with spatial structure. 
#' 
#' @param datadf The input data matrix for which GWPCA needs to be performed. It should contain numerical values only. Rows represent cells, and columns represent bioclimatic variables.
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
#' @importFrom stats princomp
#' @importFrom rio import
#' 
#' @importFrom future %<-% multisession
#' @useDynLib heterogen, .registration=TRUE
#' 
#' @return A matrix of eigenvalues
#' 
#' @examples
#' \donttest{
#' path_csv <- system.file("extdata","south.csv", package="heterogen")
#' south_csv <- rio::import(path_csv)
#' 
#' # notice: south_csv object contains x,y (lot/lat coordinates)
#' # and environmental variables
#' north_het <- gwpca_df_mc(as.matrix(south_csv), parallel = TRUE, 
#' ncores = 2, bandwidth = 0.1, tolerance = 10)
#' 
#' }
#' 
#' @export
#' 

gwpca_df_mc <- function(datadf, bandwidth = 0.2, tolerance = 5, nprocess = 10000,  parallel = FALSE, 
                     ncores = 2 , normalized = FALSE, method = 'iter', dirds='rds'){
  
  METHODS <- c("iter", "local")
  method <- match.arg(method, METHODS)
  
  message("Boosting... \n")
  
  if(class(datadf)[1] != "matrix"){
    stop('datadf, envs need to be matrix class')
  }
  
  message(paste0('Spatial setting: \n  based on bandwidth ', round(bandwidth,2)))
  
  #
  s_sample <- round(dim(datadf)[1] * (100 - tolerance) / 100, 0)
  data_s <- datadf[sample(1:dim(datadf)[1], s_sample), ]
  
  # pre files
  xy_raw <- as.matrix(datadf[,1:2])
  bg_xy <- as.matrix(data_s[,1:2])
  envs <- as.matrix(data_s[,3:dim(data_s)[2]])
  envs_t <- t(envs)
  
  # pre env
  if(normalized==FALSE){
    envin <- envs
    envin_t <- t(envin)
  }else{
    envin <- sapply(1:dim(envs)[2], function(x) (envs[, x] - min(envs[, x])) / 
                      (max(envs[, x]) - min(envs[, x])))
    
    envin_t <- t(envin)
  }
  
  if(method == 'local'){
    if(!dir.exists(dirds)) dir.create(dirds)
  }
  
  if(!parallel){
    
    steps <- seq(1, dim(xy_raw)[1], nprocess)
    kkk <- c(steps,  dim(xy_raw)[1] + 1)
    long_k <- length(kkk)
    
    heterog_env <- new.env()
    
    pasos <- 1:(length(kkk) - 1)
    pasosChar <- paste0(pasos)
    
    for (paso in pasosChar) {
      # paso <- 1
      x <- as.numeric(paso)
      heterog_env[[paso]] %<-% {
        seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
        
        in_z <- lapply(seq_rdist, function(x){
          z_in <- heterogen::gwpca_core(xy = bg_xy, p_xy = xy_raw[x, ], env = envin, env_trans = envin_t, tau = bandwidth)
          return(z_in)
        })
        
        heterog <- do.call(rbind, in_z)
        
        # save rds
        if(method == 'local'){
          saveRDS(heterog, paste0(dirds,'/het_',paso,'.rds'))
        } else{
          return(heterog)
        }
        
      }
      avance <- (x / long_k) * 100
      cat("Computation progress: ", round(avance, 1) ,"%" ,"\n")
    }
    
  if(method == 'local'){
    heterog_list <- list()
    for(i in pasosChar){
      heterog_list[[i]] <- readRDS(paste0(dirds,'/het_',i,'.rds'))
    }
    sd_scores <- do.call(rbind, heterog_list)
  }else{
    heterog_list <- as.list(heterog_env)
    heterog_names <- sort(as.numeric(as.character(names(heterog_list))))
    heterog_names <- as.character(heterog_names)
    sd_scores <- do.call(rbind, heterog_list[heterog_names])
  }
    
  }else{
    ncaval <- future::availableCores()
    
    if(ncores > ncaval){
      stop('ncores defined by user is large than available, ncores: ', ncaval)
    }
 
    # core function
    compute_het <- function(x){
      het <- heterogen::gwpca_core(xy = bg_xy, p_xy = xy_raw[x, ], env = envin, env_trans = envin_t, tau = bandwidth)
      return(het)
    }
    
    if(dim(datadf)[1] < 10000){
      
      sysinf <- Sys.info()
      os <- sysinf['sysname']
      
      if(os == "Windows") {
        cl <- parallel::makeCluster(ncores)
        tmp_zs <- parallel::parLapply(cl, 1:dim(datadf)[1], compute_het)
        
      }
      if(os == "Darwin"){
        tmp_zs <- parallel::mclapply(1:dim(datadf)[1], compute_het, mc.cores = ncores)
        
      }

      sd_scores <- do.call(rbind, tmp_zs)
      
    }else{
      sysinf <- Sys.info()
      os <- sysinf['sysname']
      
      steps <- seq(1, dim(xy_raw)[1], nprocess)
      kkk <- c(steps,  dim(xy_raw)[1] + 1)
      long_k <- length(kkk)
      
      heterog_env <- new.env()
      
      pasos <- 1:(length(kkk) - 1)
      pasosChar <- paste0(pasos)
      
      for (paso in pasosChar) {
        x <- as.numeric(paso)
        heterog_env[[paso]] %<-% {
          seq_rdist <- kkk[x]:(kkk[x + 1] - 1)
          
          if(os == "Windows") {
            cl <- parallel::makeCluster(ncores)
            z_in <- parallel::parLapply(cl, seq_rdist, compute_het)
          }
          if(os == "Darwin"){
            z_in <-  parallel::mclapply(seq_rdist, compute_het, mc.cores = ncores)
          }
          heterog <- do.call(rbind, z_in)
          
          if(method == 'local'){
            saveRDS(heterog, paste0(dirds,'/het_',paso,'.rds'))
          } 
          if (method == 'iter'){
            return(heterog)
          }
          
        }
        avance <- (x / long_k) * 100
        cat("Computation progress: ", round(avance, 1) ,"%" ,"\n")
      }
      
      if(method == 'local'){
        heterog_list <- list()
        for(i in pasosChar){
          heterog_list[[i]] <- readRDS(paste0(dirds,'/het_',i,'.rds'))
        }
        
        sd_scores <- do.call(rbind, heterog_list)
      } 
      if(method == 'iter'){
        heterog_list <- as.list(heterog_env)
        heterog_names <- sort(as.numeric(as.character(names(heterog_list))))
        heterog_names <- as.character(heterog_names)
        sd_scores <- do.call(rbind, heterog_list[heterog_names])
      }
      
    }
    
  }
  
  return(sd_scores)
  
}
