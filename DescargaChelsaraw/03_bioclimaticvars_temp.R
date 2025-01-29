#-------------------------------------------------------------------------------
# Script to generate bioclimatic variables related to temperature
# time period: 1920-1980
# Uses the tmin and tmax averages from 1920 to 1980
# BIO1, BIO2, BIO3, BIO4, BIO5, BIO6, BIO7, BIO10, BIO11.
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(furrr)
rm(list = ls())
library(magrittr)
library(stringr)
gc()
tmn_path <- list.files("tmin",recursive = T,
                       full.names = T,pattern = ".tif$")

tmx_path <- list.files("tmax",recursive = T,
                       full.names = T,pattern = ".tif$")
nona_cells <- rio::import("nona_ids.rds")

tmn_s <- str_split(tmn_path,"_",simplify = T)
n_tmn <- data.frame(year=as.numeric(tmn_s[,4]),month= as.numeric(tmn_s[,3]),
                     path=tmn_path)
n_tmn <- n_tmn[order(n_tmn$year,n_tmn$month),]
n_tmnL <- n_tmn %>% split(.$year)

tmx_s <- str_split(tmx_path,"_",simplify = T)
n_tmx <- data.frame(year=as.numeric(tmx_s[,4]),month= as.numeric(tmx_s[,3]),
                    path=tmx_path)
n_tmx <- n_tmx[order(n_tmx$year,n_tmx$month),]
n_tmxL <- n_tmx %>% split(.$year)
z=1
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

bio_temp <- seq_along(n_tmxL) %>% purrr::map(function(z){
  year_path <- file.path("bioclimatic",names(n_tmxL[z]))
  tmx_path <- n_tmxL[[z]]$path
  tmn_path <- n_tmnL[[z]]$path
  r1 <- raster::raster(tmx_path[1])
  # ------------------------------------------------------------------------------
  # mean, MIN and MAX related variables 
  # ------------------------------------------------------------------------------
  
  plan(multisession(workers = 4))
  options(future.globals.maxSize = 8500 * 1024^2)
  t_avg <- seq_along(tmn_path) %>% furrr::future_map_dfc(function(x){
    rmn <- raster::raster(tmn_path[x])
    rmn_vals <- rmn[]
    rmx <- raster::raster(tmx_path[x])
    rmx_vals <- rmx[]
    rmn_vals <- rmn_vals[nona_cells]
    rmx_vals <- rmx_vals[nona_cells]
    t_avg <-  (rmn_vals + rmx_vals)/2
    df1 <- data.frame(t_avg)
    names(df1) <- paste0("t_avg_",x)
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  t_avg <- as.matrix(t_avg)
  gc()
  bio_1_mat <- matrixStats::rowMeans2(t_avg)
  bio_1_vals <- rep(NA,raster::ncell(r1))
  bio_1_vals[nona_cells] <- bio_1_mat
  bio_1 <- r1
  bio_1[] <- round(bio_1_vals) 
  terra::writeRaster(bio_1,file.path(year_path,"bio_01.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio_1","bio_1_vals","bio_1_mat"))
  
  gc()
  # BIO 4 Temperature Seasonality (standard deviation) 
  
  bio_4_mat <- matrixStats::rowSds(t_avg) 
  bio_4_mat <- bio_4_mat*100
  bio_4_vals <- rep(NA,raster::ncell(r1))
  bio_4_vals[nona_cells] <- round(bio_4_mat)
  bio_4 <- r1
  bio_4[] <- bio_4_vals 
  terra::writeRaster(bio_4,file.path(year_path,"bio_04.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio_4","bio_4_vals","bio_4_mat"))
  gc()
  
  x <- 1:12
  xmax <- max(x)
  cuartos <- lapply(seq_along(x), function(i){
    if(i==1) return(x[1:3])
    if(i==(xmax-1)) return(c(x[i],x[i+1],x[1]))
    if(i==xmax) return(c(x[i],x[1],x[2]))
    else return(c(x[i],x[i+1],x[i+2]))
  })
  
  
  tmp_por_cuarto <- seq_along(cuartos) %>% purrr::map_dfc(function(x){
    r00 <- matrixStats::rowMeans2(t_avg[,cuartos[[x]]])
    tmp_cuarto <- data.frame(r00)
    names(tmp_cuarto) <- paste0(cuartos[[x]],collapse = "_")
    print(x)
    return(tmp_cuarto)
  })
  rm(list = c("t_avg"))
  gc()
  tmp_por_cuarto <- as.matrix(tmp_por_cuarto)
  
  ids_max <- Rfast::rowMaxs(tmp_por_cuarto)
  ids_min <- Rfast::rowMins(tmp_por_cuarto)
  tmp_cuartos_id <- data.frame(ids_max,ids_min)
  tmp_cuartos_id <- as.matrix(tmp_cuartos_id)
  rio::export(tmp_cuartos_id,file.path(year_path,"tmp_cuartos.rds"))
  
  rm(list = c("tmp_cuartos_id","ids_max","ids_min"))
  gc()
  # BIO10 Mean Temperature of Warmest Quarter 
  
  bio10_v <- rep(NA,raster::ncell(r1))
  #nona_cells <- prc_cuartos_ids[,3]
  bio10_vals <- matrixStats::rowMaxs(tmp_por_cuarto) 
  bio10_v[nona_cells] <- round(bio10_vals,2) 
  bio10 <- r1
  bio10[] <- bio10_v
  terra::writeRaster(bio10,file.path(year_path,"bio_10.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list=c("bio10","bio10_v","bio10_vals"))
  gc()
  
  bio11_v <- rep(NA,raster::ncell(r1))
  #nona_cells <- prc_cuartos_ids[,3]
  bio11_vals <- matrixStats::rowMins(tmp_por_cuarto) 
  bio11_v[nona_cells] <- round(bio11_vals) 
  bio11 <- r1
  bio11[] <- bio11_v
  terra::writeRaster(bio11,file.path(year_path,"bio_11.tif"), 
                     options = tifoptions,overwrite=TRUE)
  
  rm(list = c("bio11_v","bio11_vals","bio11"))
  gc()
  
  rm(list = c("tmp_por_cuarto"))
  gc()
  plan(multisession(workers = 4))
  options(future.globals.maxSize = 8500 * 1024^2)
  t_range <- seq_along(tmn_path) %>% furrr::future_map_dfc(function(x){
    rmn <- raster::raster(tmn_path[x])
    rmn_vals <- rmn[]
    rmx <- raster::raster(tmx_path[x])
    rmx_vals <- rmx[]
    rmn_vals <- rmn_vals[nona_cells]
    rmx_vals <- rmx_vals[nona_cells]
    t_diff <-   rmx_vals- rmn_vals
    df1 <- data.frame(t_diff)
    names(df1) <- paste0("t_diff_",x)
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  gc()
  t_range <- as.matrix(t_range)
  
  # BIO2. Mean Diurnal Range(Mean(period max-min)) 
  
  bio_2_mat <- matrixStats::rowMeans2(t_range) 
  bio_2_vals <- rep(NA,raster::ncell(r1))
  bio_2_vals[nona_cells] <- round(bio_2_mat)
  bio_2 <- r1
  bio_2[] <- bio_2_vals 
  terra::writeRaster(bio_2,file.path(year_path,"bio_02.tif"), 
                     options = tifoptions,overwrite=TRUE)
  
  rm(list = c("bio_2","bio_2_vals","bio_2_mat","t_range"))
  gc()
  
  plan(multisession(workers = 4))
  options(future.globals.maxSize = 8500 * 1024^2)
  t_max <- seq_along(tmx_path) %>% furrr::future_map_dfc(function(x){
    rmx <- raster::raster(tmx_path[x])
    rmx_vals <- rmx[]
    rmx_vals <- rmx_vals[nona_cells]
    df1 <- data.frame(rmx_vals)
    names(df1) <- paste0("t_max_",x)
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  t_max <- as.matrix(t_max)
  gc()
  # BIO5. Max Temperature of Warmest Period 
  
  bio_5_mat <- matrixStats::rowMaxs(t_max) 
  bio_5_vals <- rep(NA,raster::ncell(r1))
  bio_5_vals[nona_cells] <- round(bio_5_mat)
  bio_5 <- r1
  bio_5[] <- bio_5_vals 
  terra::writeRaster(bio_5,file.path(year_path,"bio_05.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio_5","bio_5_vals","bio_5_mat","t_max"))
  gc()
  
  plan(multisession(workers = 4))
  options(future.globals.maxSize = 8500 * 1024^2)
  t_min <- seq_along(tmn_path) %>% furrr::future_map_dfc(function(x){
    rmn <- raster::raster(tmn_path[x])
    rmn_vals <- rmn[]
    rmn_vals <- rmn_vals[nona_cells]
    df1 <- data.frame(rmn_vals)
    names(df1) <- paste0("t_min_",x)
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  t_min <- as.matrix(t_min)
  gc()
  # BIO6. Min Temperature of Coldest Period 
  
  bio_6_mat <- matrixStats::rowMins(t_min) 
  bio_6_vals <- rep(NA,raster::ncell(r1))
  bio_6_vals[nona_cells] <- round(bio_6_mat)
  bio_6 <- r1
  bio_6[] <- bio_6_vals 
  terra::writeRaster(bio_6,file.path(year_path,"bio_06.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio_6","bio_6_vals","bio_6_mat","t_min"))
  gc()
  
  bio5 <- raster::raster(file.path(year_path,"bio_05.tif"))
  bio5_vals <- bio5[]
  bio5_vals <- bio5_vals[nona_cells]
  bio6 <- raster::raster(file.path(year_path,"bio_06.tif"))
  bio6_vals <- bio6[]
  bio6_vals <- bio6_vals[nona_cells]
  bio7_vals <- bio5_vals - bio6_vals
  bio7_all <- bio5[]
  bio7_all[nona_cells] <- round(bio7_vals)
  bio7 <- bio5
  bio7[] <- bio7_all
  terra::writeRaster(bio7,file.path(year_path,"bio_07.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio5","bio5_vals","bio6","bio6_vals","bio7","bio7_vals"))
  gc()
  
  
  # BIO3. Isothermality (100 * BIO2 / BIO7) 
  bio2 <- raster::raster(file.path(year_path,"bio_02.tif"))
  bio7 <- raster::raster(file.path(year_path,"bio_07.tif"))
  #nona_cells <- rio::import("nona_ids.rds")
  bio2_v <- bio2[]
  bio2_vals <- bio2_v[nona_cells]
  bio7_v <- bio7[]
  bio7_vals <- bio7_v[nona_cells]
  bio3_v <- bio2_v  
  bio3_vals <- 100 * bio2_vals/bio7_vals
  bio3_v[nona_cells] <- round(bio3_vals,2)
  bio3 <- bio2
  bio3[] <- bio3_v
  terra::writeRaster(bio3,file.path(year_path,"bio_03.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio2","bio7","bio3","bio2_v","bio3_v","bio7_v",
              "bio2_vals","bio7_vals","bio3_vals"))
  gc()
  
})
  