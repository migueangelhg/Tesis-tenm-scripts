#-------------------------------------------------------------------------------
# Script to generate bioclimatic variables for precipitation
# time period: 2010-2016
# Uses precipitation averages from 2010 to 2016
# BIO12,BIO13,BIO14,BIO15,BIO16,BIO17
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(furrr)
rm(list = ls())
gc()
library(magrittr)

prc_path <- list.files("prec_avg_2010_2016",
                       full.names = T,pattern = ".tif$")


r1 <- raster::raster(prc_path[1])
nona_cells <- rio::import("nona_ids_2010_2016.rds")

# ------------------------------------------------------------------------------
# mean, MIN and MAX related variables 
# ------------------------------------------------------------------------------

plan(multisession(workers = 2))
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
bio12_v[nona_cells] <- bio12_vals 
bio12 <- r1
bio12[] <- bio12_v
raster::plot(bio12)
raster::writeRaster(bio12,"bios_2010_2016/bio_12.tif",overwrite=TRUE)
gc()
# BIO13. Precipitation of Wettest Period 
bio13_vals <- matrixStats::rowMaxs(prc)
bio13_v <- rep(NA,raster::ncell(r1))
bio13_v[nona_cells] <- bio13_vals 
bio13 <- r1
bio13[] <- bio13_v
raster::plot(bio13)
raster::writeRaster(bio13,"bios_2010_2016/bio_13.tif",overwrite=TRUE)
gc()
# BIO14. Precipitation of Driest Period 
bio14_vals <- matrixStats::rowMins(prc)
bio14_v <- rep(NA,raster::ncell(r1))
bio14_v[nona_cells] <- bio14_vals 
bio14 <- r1
bio14[] <- bio14_v
raster::plot(bio14)
raster::writeRaster(bio14,"bios_2010_2016/bio_14.tif",overwrite=TRUE)


#rm(list = c("r1","r1_vals","r1_ncells","nona_cells","na_cells"))
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
bio16_vals <- matrixStats::rowMins(prc_por_cuarto)
bio16_v <- rep(NA,raster::ncell(r1))
bio16_v[nona_cells] <- bio16_vals 
bio16 <- r1
raster::plot(bio16)
raster::writeRaster(bio16,"bios_2010_2016/bio_16.tif",overwrite=TRUE)
rm(list = c("bio16_vals","bio16",
            "bio16_v","bio16","bio12",
            "bio12_v,bio12_vals","bio13_v","bio13_vals",
            "bio12","bio14_v","bio14_vals"))
gc()



# P17. Precipitation of Driest Quarter 
bio17_vals <- matrixStats::rowMins(prc_por_cuarto)
bio17_v <- rep(NA,raster::ncell(r1))
bio17_v[nona_cells] <- bio17_vals 
bio17 <- r1
bio17[] <- bio17_v
raster::plot(bio17)
raster::writeRaster(bio17,"bios_2010_2016/bio_17.tif",overwrite=TRUE)
rm(list = c("bio17_vals","bio17_v","bio17"))
gc()



ids_max <- Rfast::rowMaxs(prc_por_cuarto)
ids_min <- Rfast::rowMins(prc_por_cuarto)
prc_cuartos_id <- data.frame(ids_max,ids_min)
prc_cuartos_id <- as.matrix(prc_cuartos_id)
dim(prc_cuartos_id )
rio::export(prc_cuartos_id,"prc_cuartos_ids_2010_2016.rds")
gc()

# BIO15. Precipitation Seasonality(Coefficient of Variation) 
# the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
prc <- prc+1
gc()
mean_prc <- matrixStats::rowMeans2(prc)
sd_prc <- matrixStats::rowSds(prc)
bio15_vals <- sd_prc/mean_prc
bio15_v <- rep(NA,raster::ncell(r1))
bio15_v[nona_cells] <- bio15_vals 
bio15 <- r1
bio15[] <- bio15_v
raster::plot(bio15)
gc()
raster::writeRaster(bio15,"bios_2010_2016/bio_15.tif",overwrite=TRUE)
rm(list=ls())
gc()
