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
  
  rm(list = c("r1","cell_ids_year"))
  
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