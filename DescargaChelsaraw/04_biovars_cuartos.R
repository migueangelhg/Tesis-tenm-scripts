#-------------------------------------------------------------------------------
# Script to generate bioclimatic variables for temperatures and precipitations
# of the wettest and hottest quarters
# time period: 1920-1980
# Uses the tmin and tmax averages from 1920 to 1980
# BIO8,BIO9,BIO18,BIO19
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(furrr)
library(stringr)
rm(list = ls())
gc()
library(magrittr)
library(raster)

prc_paths <- list.files("prec",
                        full.names = T,
                        pattern = ".tif$",
                        recursive = T)

prec_s <- str_split(prc_paths,"_",simplify = T)

n_prec <- data.frame(year=as.numeric(prec_s[,4]),month= as.numeric(prec_s[,3]),
                     path=prc_paths)

n_prec <- n_prec[order(n_prec$year,n_prec$month),]

n_precL <- n_prec %>% split(.$year)

tmn_path <- list.files("tmin",recursive = T,
                       full.names = T,pattern = ".tif$")

tmx_path <- list.files("tmax",recursive = T,
                       full.names = T,pattern = ".tif$")


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
nona_cells <- rio::import("nona_ids.rds")
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

secuencia <- seq_along(n_tmxL)
bio_cuartos <- secuencia[secuencia>99] %>% purrr::map(function(z){
  year_path <- file.path("bioclimatic",names(n_tmxL[z]))
  tmx_path <- n_tmxL[[z]]$path
  tmn_path <- n_tmnL[[z]]$path
  prc_path <- n_precL[[z]]$path
  prc_cuartos_ids <- rio::import(file.path(year_path,"prc_cuartos_ids.rds"))
  prc_cuartos_ids <- as.matrix(prc_cuartos_ids)
  r1 <- raster::raster(tmx_path[1])
  
  gc()
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
  
  rm(list="t_avg")
  gc()
  
  colnames(prc_cuartos_ids)
  bio8_vals <- tmp_por_cuarto[cbind(1:nrow(tmp_por_cuarto),prc_cuartos_ids[,1])]
  #rm(t_avg)
  # mat[x+nrow(mat)*(y-1)]
  bio8_v <- rep(NA,raster::ncell(r1))
  bio8_v[nona_cells] <- bio8_vals 
  bio8 <- r1
  bio8[] <- bio8_v
  terra::writeRaster(bio8,file.path(year_path,"bio_08.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio8","bio8_vals","bio8_v"))
  gc()
  
  # BIO9. Mean Temperature of Driest Quarter 
  #t_avg <- as.matrix(t_avg)
  #colnames(prc_cuartos_ids)
  bio9_vals <- tmp_por_cuarto[cbind(1:nrow(tmp_por_cuarto),
                                    prc_cuartos_ids[,2])]
  bio9_v <- rep(NA,raster::ncell(r1))
  bio9_v[nona_cells] <- bio9_vals 
  bio9 <- r1
  bio9[] <- bio9_v
  terra::writeRaster(bio9,file.path(year_path,"bio_09.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio9","bio9_vals",
              "bio9_v","prc_cuartos_ids",
              "tmp_por_cuarto"))
  gc()
  
  
  tmp_cuartos_ids  <- rio::import(file.path(year_path,"tmp_cuartos.rds"))
  
  plan(multisession(workers = 4))
  options(future.globals.maxSize = 8500 * 1024^2)
  prc <- seq_along(prc_path) %>% furrr::future_map_dfc(function(x){
    pre <- raster::raster(prc_path[x])
    pre_vals <- pre[]
    pre_vals <- pre_vals[nona_cells]
    df1 <- data.frame(pre_vals)
    names(df1) <- paste0("pre",x)
    return(df1)
  },.progress = TRUE)
  plan(sequential)
  prc <- as.matrix(prc)
  gc()
  
  x <- 1:12
  xmax <- max(x)
  
  cuartos <- lapply(seq_along(x), function(i){
    if(i==1) return(x[1:3])
    if(i==(xmax-1)) return(c(x[i],x[i+1],x[1]))
    if(i==xmax) return(c(x[i],x[1],x[2]))
    else return(c(x[i],x[i+1],x[i+2]))
  })
  
  
  
  prc_por_cuarto <- seq_along(cuartos) %>% purrr::map_dfc(function(x){
    r00 <- matrixStats::rowSums2(prc[,cuartos[[x]]])
    prc_cuarto <- data.frame(r00)
    names(prc_cuarto) <- paste0(cuartos[[x]],collapse = "_")
    print(x)
    return(prc_cuarto)
  })
  rm("prc")
  gc()
  
  # BIO18. Precipitation of Warmest Quarter 
  
  bio18_vals <- prc_por_cuarto[cbind(1:nrow(tmp_cuartos_ids),
                                     tmp_cuartos_ids[,1])]
  #rm(t_avg)
  # mat[x+nrow(mat)*(y-1)]
  bio18_v <- rep(NA,raster::ncell(r1))
  bio18_v[nona_cells] <- bio18_vals 
  bio18 <- r1
  bio18[] <- bio18_v
  #rm(list = c("prc_cuartos_ids","bio8_vals","bio8_v"))
  gc()
  terra::writeRaster(bio18,file.path(year_path,"bio_18.tif"), 
                     options = tifoptions,overwrite=TRUE)
  
  
  bio19_vals <- prc_por_cuarto[cbind(1:nrow(tmp_cuartos_ids),
                                     tmp_cuartos_ids[,2])]
  gc()
  #rm(t_avg)
  # mat[x+nrow(mat)*(y-1)]
  bio19_v <- rep(NA,raster::ncell(r1))
  bio19_v[nona_cells] <- bio19_vals 
  bio19 <- r1
  bio19[] <- bio19_v
  gc()
  terra::writeRaster(bio19,file.path(year_path,"bio_19.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio18","bio18_vals","bio18_v",
              "bio19","bio19_vals","bio19_v",
              "prc_por_cuarto",
              "tmp_cuartos_ids"))
  print(year_path)
  gc()
  return()
})
