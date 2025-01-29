#-------------------------------------------------------------------------------
# Script to generate bioclimatic variables for precipitation
# time period: 1920-1980
# Uses precipitation averages from 1920 to 1980
# BIO12,BIO13,BIO14,BIO15,BIO16,BIO17
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(furrr)
library(stringr)
rm(list = ls())
gc()
library(magrittr)

prc_paths <- list.files("prec",
                        full.names = T,
                        pattern = ".tif$",
                        recursive = T)

prec_s <- str_split(prc_paths,"_",simplify = T)

n_prec <- data.frame(year=as.numeric(prec_s[,4]),month= as.numeric(prec_s[,3]),
                     path=prc_paths)

n_prec <- n_prec[order(n_prec$year,n_prec$month),]

n_precL <- n_prec %>% split(.$year)

nona_cells <- rio::import("nona_ids.rds")
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")


bio_prec <- seq_along(n_precL) %>% purrr::map(function(z){
  year_path <- file.path("bioclimatic",names(n_precL[z]))
  prc_path <- n_precL[[z]]$path
  r1 <- raster::raster(prc_path[1])
  # ------------------------------------------------------------------------------
  # mean, MIN and MAX related variables 
  # ------------------------------------------------------------------------------
  
  plan(multisession(workers = 8))
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
  # BIO12. Annual Precipitation 
  bio12_vals <- matrixStats::rowSums2(prc) 
  bio12_v <- rep(NA,raster::ncell(r1))
  bio12_v[nona_cells] <- as.integer(round(bio12_vals)) 
  bio12 <- r1
  bio12[] <- bio12_v
  #bio12 <- as.integer(bio12)
  
  terra::writeRaster(bio12,file.path(year_path,"bio_12.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio12","bio12_vals","bio12_v"))
  gc()
  # BIO13. Precipitation of Wettest Period 
  bio13_vals <- matrixStats::rowMaxs(prc)
  bio13_v <- rep(NA,raster::ncell(r1))
  bio13_v[nona_cells] <- as.integer(round(bio13_vals)) 
  bio13 <- r1
  bio13[] <- bio13_v
  terra::writeRaster(bio13,file.path(year_path,"bio_13.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio13","bio13_vals","bio13_v"))
  gc()
  
  # BIO14. Precipitation of Driest Period 
  
  bio14_vals <- matrixStats::rowMins(prc) 
  bio14_v <- rep(NA,raster::ncell(r1))
  bio14_v[nona_cells] <- as.integer(round(bio14_vals))
  bio14 <- r1
  bio14[] <- bio14_v
  
  
  terra::writeRaster(bio14,file.path(year_path,"bio_14.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio14","bio14_vals","bio14_v"))
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
  #rm(list = c("prc"))
  
  prc_por_cuarto <- as.matrix(prc_por_cuarto)
  gc()
  
  # P16. Precipitation of Wettest Quarter 
  bio16_vals <- matrixStats::rowMaxs(prc_por_cuarto) 
  bio16_v <- rep(NA,raster::ncell(r1))
  bio16_v[nona_cells] <- as.integer(round(bio16_vals)) 
  bio16 <- r1
  bio16[] <- bio16_v
  
  terra::writeRaster(bio16,file.path(year_path,"bio_16.tif"), 
                     options = tifoptions,overwrite=TRUE)
  
  rm(list = c("bio16_vals","bio16",
              "bio16_v","bio16"))
  gc()
  
  # P17. Precipitation of Driest Quarter 
  bio17_vals <- matrixStats::rowMins(prc_por_cuarto)
  bio17_v <- rep(NA,raster::ncell(r1))
  bio17_v[nona_cells] <- as.integer(round(bio17_vals)) 
  bio17 <- r1
  bio17[] <- bio17_v
  
  terra::writeRaster(bio17,file.path(year_path,"bio_17.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = c("bio17_vals","bio17_v","bio17"))
  gc()
  
  ids_max <- Rfast::rowMaxs(prc_por_cuarto)
  ids_min <- Rfast::rowMins(prc_por_cuarto)
  prc_cuartos_id <- data.frame(ids_max,ids_min)
  prc_cuartos_id <- as.matrix(prc_cuartos_id)
  dim(prc_cuartos_id )
  rio::export(prc_cuartos_id,file.path(year_path,"prc_cuartos_ids.rds"))
  rm(list = "prc_por_cuarto")
  gc()
  
  prc <- prc+1
  gc()
  mean_prc <- matrixStats::rowMeans2(prc)
  sd_prc <- matrixStats::rowSds(prc)
  bio15_vals <- round(sd_prc/mean_prc,2)
  bio15_v <- rep(NA,raster::ncell(r1))
  bio15_v[nona_cells] <- bio15_vals 
  bio15 <- r1
  bio15[] <- bio15_v
  raster::plot(bio15)
  #gc()
  terra::writeRaster(bio15,file.path(year_path,"bio_15.tif"), 
                     options = tifoptions,overwrite=TRUE)
  rm(list = ls())
  gc()
  
  return()
  
})

