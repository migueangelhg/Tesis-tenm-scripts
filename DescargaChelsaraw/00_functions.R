ex_by_year2 <- function(this_species,layers_pattern){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  years <- sort(unique(this_species$sp_occs_year))
  capasByRes <- lapply(this_species$layers_path_by_year,
                       function(x){
                         x[grep(layers_pattern,x = x)]
                       })
  
  years_env_fun <- function(x){
    year_env <- x$years_in_occs[1]
    env_layers <- raster::stack(capasByRes[[paste0(year_env)]])
    env_sp <- raster::extract(x = env_layers,
                              y = x[,this_species$lon_lat_vars],
                              df=TRUE)
    
    names(env_sp) <- c("ID_YEAR",
                       names(env_sp)[2:length(names(env_sp))])
    coords.x1 <- x[,this_species$lon_lat_vars[1]]
    coords.x2 <- x[,this_species$lon_lat_vars[2]]
    
    env_sp$ID_YEAR <- year_env
    n_70 <- floor(dim(env_sp)[1]*.7)
    env_sp$trian_test <- NA
    train_index <- sample(1:dim(env_sp)[1],n_70)
    env_sp$trian_test[train_index] <- "Train"
    env_sp$trian_test[-train_index] <- "Test"
    return(cbind(coords.x1,coords.x2,env_sp))
  }
  
  df_occs_year <- this_species$oocs_data
  df_occs_yearL <- df_occs_year %>%
    split(.$years_in_occs)
  nyearss <- length(df_occs_yearL)
  pb <- utils::txtProgressBar(min = 0, max = nyearss, initial = 0,style = 3) 
  years_env <- seq_along(df_occs_yearL ) %>%
    purrr::map_df(function(x){
      setTxtProgressBar(pb,x)
      res <- years_env_fun(df_occs_yearL[[x]])
      return(res)
    })

  names(years_env)[1:2] <- this_species$lon_lat_vars
  
  train_data <- which(years_env$trian_test=="Train")
  sp.temp.data.env <- list(sp_coords = years_env[,this_species$lon_lat_vars],
                           coords_env_data_all = years_env,
                           env_data_train = years_env[train_data,
                                                      4:(dim(years_env)[2]-1)],
                           env_data_test = years_env[-train_data,
                                                     4:(dim(years_env)[2]-1)],
                           test_data =years_env[-train_data,],
                           sp_occs_year = years_env$ID_YEAR,
                           oocs_data = this_species$oocs_data,
                           lon_lat_vars = this_species$lon_lat_vars,
                           layers_path_by_year = capasByRes)
  class(sp.temp.data.env) <- c("list", "sp.temporal.modeling","sp.temporal.env")
  
  
  return(sp.temp.data.env)
  
}


library(furrr)
ex_by_year <- function(this_species,layers_pattern,ncores=10,train_prop=0.7){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  years <- sort(unique(this_species$sp_occs_year))
  capasByRes <- lapply(this_species$layers_path_by_year,
                       function(x){
                         x[grep(layers_pattern,x = x)]
                       }) 
  
  capasByResDF <- seq_along(capasByRes) %>% purrr::map_df(function(x){
    tdf <- data.frame(time=names(capasByRes[x]),
                      layer_path=capasByRes[[x]])
  })
  
  cell_ids_year <- raster::cellFromXY(raster::raster(capasByResDF$layer_path[1]),
                                      this_species$oocs_data[,this_species$lon_lat_vars])
  
  dfy <- data.frame(cell_ids_year,ID_YEAR = this_species$oocs_data$year,
                    this_species$oocs_data[,this_species$lon_lat_vars])
  plan(tweak("multisession",workers=ncores))
  ex_time <- 1:nrow(capasByResDF) %>% furrr::future_map_dfr(function(x){
    time_obs <- dfy  %>% dplyr::filter(ID_YEAR == capasByResDF$time[x])
    env_layers <- raster::raster(capasByResDF$layer_path[x])
    layer_val <- env_layers[time_obs$cell_ids_year]
    df1 <- data.frame(time_obs[,c(3,4,2)],
                      #ID_YEAR = capasByResDF$time[x],
                      layer_val,var_name = names(env_layers))
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  gc()
  years_env <- tidyr::pivot_wider(ex_time, 
                                  names_from = var_name, 
                                  values_from = layer_val) 
  years_envL <- years_env %>% split(.$ID_YEAR)
  
  trian_test <- seq_along(years_envL) %>% purrr::map(function(x){
    ndata <- nrow(years_envL[[x]]) 
    if(ndata==1) train_test <- "Train"
    if(ndata==2) train_test <- c("Train","Test")
    if(ndata==3) {
      train_test <- sample(c("Train","Train","Test"),size = 3)
    } else{
      train_test <- rep("Test",ndata)
      ids_train <- sample(ndata,size = ceiling(ndata*train_prop))
      train_test[ids_train] <- "Train"
      
    }
    return(train_test)
  }) %>% unlist()
  years_env$trian_test <- trian_test
  
  train_data <- which(years_env$trian_test=="Train")
  sp.temp.data.env <- list(sp_coords = years_env[,this_species$lon_lat_vars],
                           coords_env_data_all = years_env,
                           env_data_train = years_env[train_data,
                                                      4:(dim(years_env)[2]-1)],
                           env_data_test = years_env[-train_data,
                                                     4:(dim(years_env)[2]-1)],
                           test_data =years_env[-train_data,],
                           sp_occs_year = years_env$ID_YEAR,
                           oocs_data = this_species$oocs_data,
                           lon_lat_vars = this_species$lon_lat_vars,
                           layers_path_by_year = capasByRes)
  class(sp.temp.data.env) <- c("list", "sp.temporal.modeling","sp.temporal.env")
  
  
  return(sp.temp.data.env)
  
}


bg_by_year <- function(this_species,layers_pattern,ncores=10,
                       buffer_ngbs=10,n_bg=50000){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  years <- sort(unique(this_species$sp_occs_year))
  capasByRes <- lapply(this_species$layers_path_by_year,
                       function(x){
                         x[grep(layers_pattern,x = x)]
                       }) 
  
  capasByResDF <- seq_along(capasByRes) %>% purrr::map_df(function(x){
    tdf <- data.frame(time=names(capasByRes[x]),
                      layer_path=capasByRes[[x]])
  })
  r1 <- raster::raster(capasByResDF$layer_path[1])
  cell_ids_year <- raster::cellFromXY(r1,
                                      this_species$oocs_data[,this_species$lon_lat_vars])
  rm(r1)
  
  dfy <- data.frame(cell_ids_year,ID_YEAR = this_species$oocs_data$year,
                    this_species$oocs_data[,this_species$lon_lat_vars])
  
  df_samps <-  dfy %>% dplyr::group_by(ID_YEAR) %>% 
    dplyr::summarise(samp_prop = n()/nrow(dfy),
                     n_samples = ceiling(samp_prop*n_bg))
  paths_layers <-   capasByResDF %>% dplyr::group_by(time) %>% dplyr::sample_n(1)
  names(paths_layers)[1] <- "ID_YEAR"
  paths_layers$ID_YEAR <- as.numeric(paths_layers$ID_YEAR)
  df_samps <- df_samps %>% dplyr::inner_join(paths_layers)
  
  nbase <- 2*buffer_ngbs+1
  ngMat <- base::matrix(rep(1,nbase*nbase),
                        ncol =nbase,byrow = T )
  ngMat[buffer_ngbs+1,buffer_ngbs+1] <- 0
  
  
  ddL <- dfy %>% split(.$ID_YEAR)
  cells_to_samp <-  seq_along(ddL)  %>% purrr::map(function(z){
    r1 <- raster::raster(df_samps$layer_path[z])
    adj_cells <- raster::adjacent(x = r1,cells=ddL[[z]]$cell_ids_year,
                                  directions = ngMat)
    rcells <- unique(adj_cells[,2])
    sam <- sample(rcells,size = df_samps$n_samples[z])
    return(sam)
  })
  names(cells_to_samp) <- names(ddL)
  
  plan(tweak("multisession",workers=ncores))
  ex_time <- 1:nrow(capasByResDF) %>% furrr::future_map_dfr(function(x){
    cellids <-  cells_to_samp[[as.character(capasByResDF$time[x])]]
    #buffers_year[[as.character(time_obs$ID_YEAR[1])]]
    env_layers <- raster::raster(capasByResDF$layer_path[x])
    layer_val <- na.omit(env_layers[cellids])
    df1 <- data.frame(ID_YEAR = capasByResDF$time[x],
                      layer_val,var_name = names(env_layers))
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  gc()
  bg_env <- tidyr::pivot_wider(ex_time, 
                               names_from = var_name, 
                               values_from = layer_val) 
  
  bg_env <-   seq_along(bg_env$ID_YEAR) %>% purrr::map_dfr(function(x){
    ID_YEAR <- rep(bg_env$ID_YEAR[[x]],length(bg_env[[2]][[x]]))
    df_year <- seq_along(bg_env[-1]) %>% purrr::map_dfc(function(y){
      data <- bg_env[-1]
      value <- data[[y]][[x]]
      df1 <- data.frame(value)
      names(df1) <- names(data[y])
      return(df1)
    })
    df_res <- data.frame(ID_YEAR,df_year)
    return(df_res)
  })
  
  return( bg_env)
  sp.temp.data.env <- list(sp_coords = years_env[,this_species$lon_lat_vars],
                           coords_env_data_all = years_env,
                           env_data_train = years_env[train_data,
                                                      4:(dim(years_env)[2]-1)],
                           env_data_test = years_env[-train_data,
                                                     4:(dim(years_env)[2]-1)],
                           test_data =years_env[-train_data,],
                           sp_occs_year = years_env$ID_YEAR,
                           oocs_data = this_species$oocs_data,
                           lon_lat_vars = this_species$lon_lat_vars,
                           layers_path_by_year = capasByRes,
                           env_bg = bg_env)
  class(sp.temp.data.env) <- c("list", "sp.temporal.modeling","sp.temporal.env")
  
  
  return(sp.temp.data.env)
  
}
#' Function to find the best n-dimensional ellipsoid model using Partial Roc as a performance criteria.
#' @param this_species, Species Temporal Environment "sp.temporal.env" object see \code{\link[hsi]{extract_by_year}}.
#' @param cor_threshold Threshold valuefrom which it is considered that the correlation is high see \code{\link[hsi]{correlation_finder}}.
#' @param nbg_points Number of background points used to compute partial ROC test. See \code{\link[ntbox]{sample_envbg}} for more details.
#' @param omr_criteria Omission rate used to select best models. See \code{\link[ntbox]{ellipsoid_selection}} for more details.
#' @param ellipsoid_level The proportion of points to be included inside the ellipsoid see \code{\link[hsi]{ellipsoidfit}}.
#' @param nvars_to_fit Number of variables that will be used to model.
#' @param E  Amount of error admissible for Partial Roc test (by default =.05). Value should range between 0 - 1. see \code{\link[hsi]{PartialROC}}
#' @param RandomPercent Occurrence points to be sampled in randomly for the boostrap of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @param NoOfIteration Number of iteration for the bootstrapping of the Partial Roc test \code{\link[hsi]{PartialROC}}.
#' @param parallel Logical argument to run computations in parallel. Default TRUE
#' @param n_cores Number of cores to be used in parallelization. Default 4
#' @return A "sp.temp.best.model" object with metadata of the best model given the performance of the Partial Roc test.
#' @export

find_best_model_ntbox2 <- function(this_species,cor_threshold=0.9,
                                  nbg_points=50000,
                                  omr_criteria =0.1,
                                  ellipsoid_level=0.975,
                                  nvars_to_fit=c(2,3),
                                  E = 0.05,
                                  RandomPercent = 50,
                                  NoOfIteration=1000,
                                  parallel=TRUE,
                                  n_cores=6){
  #ntbox_funcs <- system.file("helpers/helpers_from_ntbox.R",package = "hsi")
  #source(ntbox_funcs)
  stopifnot(inherits(this_species, "sp.temporal.env"))
  n_nas <- floor(dim(this_species$env_data_train)[1]*0.1)
  env_train <- this_species$env_data_train
  env_test <-this_species$env_data_test
  rm_layers <- unlist(sapply( 1:dim(env_train)[2], function(x){
    if(length(which(is.na(env_train[,x]))) > n_nas) return(x)
  } ))
  
  if(!is.null(rm_layers)){
    env_train <- stats::na.omit(env_train[,-rm_layers])
    env_test <- stats::na.omit(env_test[,-rm_layers])
  }
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("The total number of occurrence records that will be used for model validation is:",
      nrow(env_test), "\n")
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
  numericIDs <- which(sapply(env_train, is.numeric))
  cor_matrix <- stats::cor(na.omit(env_train[,numericIDs]))
  
  find_cor   <- correlation_finder(cor_mat = cor_matrix,
                                   threshold = cor_threshold,
                                   verbose = F)
  cor_filter <- find_cor$descriptors
  year_to_search <- min(as.numeric(names(this_species$layers_path_by_year)))
  
  env_layers <- raster::stack(this_species$layers_path_by_year[[paste0(year_to_search)]])
  bg_samp <- sample_envbg(env_layers,nbg = nbg_points)
  ids_fit <- which(nvars_to_fit>length(cor_filter))
  if(length(ids_fit)>0L)
    nvars_to_fit <- nvars_to_fit[-ids_fit]
  if(length(nvars_to_fit) ==0L)
    nvars_to_fit <- 2:length(cor_filter)
  plan(sequential)
  seleccion_mods <- hsi::ellipsoid_selection(env_train = env_train,
                                             env_test = env_test,
                                             env_vars = cor_filter,
                                             nvarstest = nvars_to_fit,
                                             level = ellipsoid_level,
                                             mve = TRUE,
                                             ncores = n_cores,
                                             comp_each = 100,
                                             env_bg = bg_samp,
                                             parallel = T,
                                             omr_criteria = omr_criteria,
                                             proc = TRUE,
                                             proc_iter = NoOfIteration,
                                             rseed = 111)
  best_mod_vars <- stringr::str_split(seleccion_mods$fitted_vars,",")
  env_all <- na.omit(this_species$coords_env_data_all[,best_mod_vars[[1]]])
  best_model_metadata <- cov_center(data = env_all,
                                    level = ellipsoid_level,
                                    mve = TRUE,
                                    vars = best_mod_vars[[1]])
  
  
  sp.temp.best.model <- list(sp_coords = this_species$sp_coords,
                             coords_env_data_all = this_species$coords_env_data_all,
                             env_data_train = this_species$env_data_train,
                             env_data_test = this_species$env_data_test,
                             test_data = this_species$test_data,
                             sp_occs_year = this_species$sp_occs_year,
                             oocs_data = this_species$oocs_data,
                             lon_lat_vars = this_species$lon_lat_vars,
                             layers_path_by_year = this_species$layers_path_by_year,
                             best_model_metadata= best_model_metadata,
                             model_selection_results =seleccion_mods,
                             ellipsoid_level =ellipsoid_level)
  class(sp.temp.best.model) <- c("list", "sp.temporal.modeling","sp.temporal.env","sp.temp.best.model")
  
  
  return(sp.temp.best.model)
  
}