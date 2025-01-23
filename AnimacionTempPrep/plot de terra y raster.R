library(raster)
library(ggplot2)
library(terra)
t <-mapa1975 <- raster::raster("/media/miguel/TOSHIBA EXT/Biomex/1986/biomex_bio_01.tif")
r <- rast(mapa1975)

### usar función extent si es posible

#e <- extent(mapa1975)

terra::plot(r, type = "interval")
ex <- c(-119,101,10,40)

plot(r, plg=list(ext=ex, title="Title\n", title.cex=1.25),
     pax=list(side=1:4, retro=F))
north("topleft")

raster::plot(t)
png::writePNG(mapa1975)



p <- raster::raster("/media/miguel/TOSHIBA EXT/Biomex/1986/biomex_bio_12.tif")

plot(t,p)
