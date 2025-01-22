library(raster)
library(rgdal)
library(dplyr)
library(sp)
library(sf)
library(tidyverse)
library(gtools)
library(rgeos)

rm(list=ls())
getwd()
setwd("E:/")




años <- list.dirs("bioclimatic/")

archivos <- list.files(años, pattern = ".tif$", full.names = TRUE) %>%
  mixedsort()

capas <- lapply(archivos, FUN= raster) %>%
  stack()
  
mex <- rgdal::readOGR("E:/mex_gadm/gadm36_MEX_0.shp")

plan(multisession)
capas_mex <- raster::crop(capas, mex) %>%
  raster::mask(mex)
capas_mex <- unstack(capas_mex)
nombresraster <- paste0("biomex",1:19)

writeRaster()



Map("writeRaster", x = capas_mex, filename = paste0("E:/capas_bioclim_mex/1901/", nombresraster, ".tif"))


bi_clim <- raster::raster("E:/capas_bioclim_mex/1901/biomex19.tif")
plot(bi_clim)
