# Las capas biclimaticas son archivos tipo .tif
#Tagged Image File Format/formato de archivo de imagen etiquetada
#formato de archivo de imagen para almacenar raster  
# Imagenes de mapa de bits ((raster))= (Imagen bidimensional con información
#una matriz rectangular)
setwd("E:/bioclimatic")
library(raster)
library(sp)
library(rgdal)
install.packages("rgeos")
library(tidyverse)
install.packages("gtools")
library(gtools)
library(rgeos)
install.packages("tidyverse")
install.packages("GDAL")
library(gdal)
normalizePath(E://bioclimatic)

rm(list=ls())



files <- list.files("1901", full.names = TRUE, pattern = ".tif$") %>%
  mixedsort()

layers <- lapply(files, FUN = raster) %>%
  stack()

mex <- shapefile("E:/mex_gadm/gadm36_MEX_0.shp")

plot(layers[[1]])
plot(mex, add = TRUE)

layers_mex <- raster::crop(layers, mex) %>%
  raster::mask(mex)
layers_mex <- unstack(layers_mex)
namesRaster <- paste0 ("biomex", 1:19)
Map("writeRaster", x = layers_mex, filename = paste0("E:/capas_bioclim_mex/1901/", namesRaster, ".tif"))



bi_clim <- rgdal::(dsn = "E:\bioclimatic\1901\bio_019.tif")

bi_clim <- raster::raster("E:/capas_bioclim_mex/1901/biomex1.tif")
plot(bi_clim)
str(bi_clim)
class(bi_clim)
?readall

par("mar") 

par(mar=c(1,1,1,1))  
getwd()
getwd()
