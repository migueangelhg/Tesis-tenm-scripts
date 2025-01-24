#-------------------------------------------------------------------------------
# Script to create tmin averages
# time period: 2010-2016
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(raster)
library(stringr)
library(matrixStats)
library(furrr)
rm(list = ls())
gc()
tmin_paths <- list.files("tmin_2010_2016",pattern = ".tif$",full.names = TRUE)
tmin_2010_2016 <- grep("_20[1-2][0-9]_V.1.0.tif",tmin_paths)

tmin_paths_n <- tmin_paths[tmin_2010_2016]
## nNA value -32768
r1 <- raster::raster(tmax_paths_n[1])
r1_vals <- rep(NA,raster::ncell(r1))
nona_cells <- rio::import("nona_ids_2010_2016.rds")

# Data 
tmin_s <- str_split(tmin_paths_n,"_",simplify = T)

n_tmin <- data.frame(year=as.numeric(tmin_s[,6]),month= as.numeric(tmin_s[,5]),
                     path=tmin_paths_n)
n_tminL <- n_tmin %>% split(.$month)
x=1
if(!dir.exists("tmin_avg_2010_2016")) dir.create("tmin_avg_2010_2016")
months <- c(paste0("0",1:9),10:12)
filenames <- paste0("tmin_",months,".tif")
filepaths <- file.path("tmin_avg_2010_2016",filenames)

r_tmin <- seq_along(n_tminL) %>% purrr::map(function(x){
  df0 <- n_tminL[[x]]
  nyears <- 1:nrow(df0)
  npasos <- 6
  paso <- round(nrow(df0)/npasos)
  pasos <- lapply(1:npasos, function(h){
    if(h==1) return(1:paso)
    if(h==npasos) return((paso*(h-1)+1):nrow(df0))
    else return((paso*(h-1)+1):(paso*h))
  })
  rrange <- seq_along(pasos) %>% purrr::map_dfc(function(y){
    nyears <- pasos[[y]]
    plan(multisession(workers = 4))
    options(future.globals.maxSize = 8500 * 1024^2)
    vals_df <- nyears %>% furrr::future_map_dfc(function(z){
      r0 <- raster::raster(df0$path[z])
      r0_vals <- r0[]
      var_name <- names(r0)
      rm(r0)
      r0_nonas <- r0_vals[nona_cells]
      rm(r0_vals)
      df_vals <- data.frame(rval = r0_nonas)
      names(df_vals) <- var_name
      return(df_vals)
    },.progress = TRUE)
    plan(sequential)
    vals_df <- as.matrix(vals_df)
    mean_val <- matrixStats::rowMeans2(vals_df)
    rdf <- data.frame(mean_val)
    return(rdf)
  })
  #vals <- matrixStats::rowSds(as.matrix(rrange))
  #sss <- which(vals >1)
  #nnn <- rrange[sss,]
  rrange <- as.matrix(rrange)
  mean_all <- matrixStats::rowMeans2(rrange)
  r1_vals[nona_cells] <- mean_all
  r1[] <- r1_vals
  raster::writeRaster(r1,filepaths[x],overwrite=TRUE)
  print(x)
  return()
})
rm(list = ls())
gc()