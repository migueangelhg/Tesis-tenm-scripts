library(furrr)
ex_by_year <- function(this_species,layers_pattern,ncores=10,train_prop=0.7,buffer_distance=1){
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
  
  #raster::adjacent(raster::raster(capasByResDF$layer_path[1]),cell_ids_year)

  buffers_year <-  this_species$oocs_data %>% split(.$year) %>% purrr::map(function(x){
    time_obs <- x
    time_obs_sp <- sp::SpatialPointsDataFrame(time_obs[,3:4],
                                              data = time_obs,
                                              proj4string = raster::crs(r1))
    
    bf_bg <- raster::buffer(time_obs_sp,width=buffer_distance*1000)
    return(bf_bg)
  })
  
  
  dfy <- data.frame(cell_ids_year,ID_YEAR = this_species$oocs_data$year,
                    this_species$oocs_data[,this_species$lon_lat_vars])
  plan(tweak("multisession",workers=ncores))
  ex_time <- 1:nrow(capasByResDF) %>% furrr::future_map_dfr(function(x){
    time_obs <- dfy  %>% dplyr::filter(ID_YEAR == capasByResDF$time[x])
    #buffers_year[[as.character(time_obs$ID_YEAR[1])]]
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