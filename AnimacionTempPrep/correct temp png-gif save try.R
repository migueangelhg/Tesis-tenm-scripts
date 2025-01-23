library(raster)
library(dplyr)
library(rasterVis)
library(future)
library(pdftools)
library(magick)


path_bio1 <- "/media/miguel/TOSHIBA EXT/Temperatura_animacion"
if(!dir.exists(path_bio1))dir.create(path_bio1)
dirs <- list.dirs("/media/miguel/TOSHIBA EXT/Biomex/")
anios <-list.files("/media/miguel/TOSHIBA EXT/Biomex/",
                   full.names = FALSE)
fils <- list.files(dirs, pattern = "bio_01.tif$", full.names = TRUE)



plan(multisession(workers = 10))
x=1

anima <- seq_along(fils)%>% furrr::future_map(function(x){
  r <- raster::raster(fils[x])
  r <- r/10
  png(paste0(path_bio1,"/",anios[x],".png"),
      width = 1200,
      height = 1000,
      res = 215)
  pngs <- levelplot(r, main = paste0("Temperatura promedio anual ",anios[x]),
                    margin=FALSE,
                    colorkey=list(space="bottom"),
                    par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
  print(pngs)
  dev.off()
})
plan(sequential)

library(gifski)
png_files <- list.files(path = path_bio1, pattern = ".*png$", full.names = TRUE)
gifski(png_files,
       gif_file = "temp_anima.gif",
       width = 1200,
       height = 1000,
       delay = .3)

