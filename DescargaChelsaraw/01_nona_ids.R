#-------------------------------------------------------------------------------
# Script to identify non-nas pixels
# Indices will be used to better manage memory
# time period: 1901-2016
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------


library(furrr)
rm(list = ls())
gc()
library(magrittr)

prc_path <- list.files("prec",
                       full.names = T,pattern = ".tif$")
tmn_path <- list.files("tmin",
                       full.names = T,pattern = ".tif$")

tmx_path <- list.files("tmax",
                       full.names = T,pattern = ".tif$")

rtmnp <- raster::raster(tmn_path[1])
rprc <- raster::raster(prc_path[1])
rtmxp <- raster::raster(tmx_path[1])

na_val <- raster::cellStats(rtmnp,min)
raster::NAvalue(rtmnp) <- na_val
raster::NAvalue(rprc) <- na_val
raster::NAvalue(rtmxp) <- na_val

rtmnp_vals <- rtmnp[]
rtmxp_vals <- rtmxp[]
rprc_vals <- rprc[]
nona_tmxp <- which(!is.na(rtmxp_vals))
nona_tmnp <- which(!is.na(rtmnp_vals))
nona_prc <- which(!is.na(rprc_vals))
nonas_1 <- intersect(nona_tmnp,nona_prc)
nonas_2 <- intersect(nonas_1,nona_tmxp)
length(nonas_1)
length(nonas_2)

nona_mat <- matrix(nonas_2,ncol = 1)
mode(nona_mat) <- "integer"
class(nona_mat[,1])

rio::export(nona_mat[,1],"nona_ids.rds")
