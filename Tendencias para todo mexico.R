library(rgdal)
library(sf)
library(raster)
library(dplyr)
library(furrr)
rm(list = ls())
gc()
anps <- readOGR("E:/poligonos_anps_mex/WDPA_WDOECM_Nov2021_Public_MEX_shp-polygons.shp")
plot(anps)


capas_path <- list.files("E:/Biomex", pattern = "biomex_bio_01.tif$",recursive = T, full.names = T)
r <- raster("E:/Biomex/biomex_2013/biomex_bio_01.tif")
plot(r)
r<- raster::stack(capas_path)
names(r) <- paste0("bio_01_", 1901:2016)
#cell_ids_anp <- unlist(raster::cellFromPolygon(r,anps))

#xy <- sp::coordinates(r[[1]])[cell_ids_anp,]

#datos <- data.frame (xy,r[cell_ids_anp])
#datos <- na.omit(datos)
#datos[1,-(1:2)]



#hola <- data.frame(Temperatura = as.numeric(datos[1,-(1:2)]),anio=1901:2016)


cellvals <- r[[1]][]
cellIDs <-  which(!is.na(cellvals))
XY <- sp::coordinates(r[[1]])[cellIDs,]
datos <- data.frame(XY, r[cellIDs])

plan(multisession(workers  = 10))
options(future.globals.maxSize = 2500 * 1024^2)
df1 <- 1:nrow(datos) %>% furrr::future_map_dfr(function(x){
  hola <- data.frame(Temperatura = as.numeric(datos[x,-(1:2)]),anio=1901:2016)
  m <- lm(Temperatura ~ anio,data = hola)
  resultado <- summary(m)
  pend <- coef(m)[2]
  Pval <- resultado$coefficients[2,4]
  temperaturas <- data.frame(datos[x,1:2],pend,Pval)
  return(temperaturas)
},.progress = T) 

plan(sequential)

r2 <- r[[1]]*NA
cellids <- cellFromXY(r2,df1[,1:2])
r2[cellids] <- df1$pend
plot(r2)
