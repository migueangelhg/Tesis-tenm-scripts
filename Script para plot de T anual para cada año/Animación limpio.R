library(raster)
library(dplyr)
library(rasterVis)
library(future)
#install.packages("magick")
#install.packages("pdftools")
library(pdftools)
library(magick)
#library(ggplot2)

####
path_anima <- "/media/miguel/TOSHIBA EXT/imagenes_para_animacion"
if(!dir.exists(path_anima))dir.create(path_anima)
####
path_bio1 <- "/media/miguel/TOSHIBA EXT/Temperatura_animacion"
if(!dir.exists(path_bio1))dir.create(path_bio1)
path_prep <- "/media/miguel/TOSHIBA EXT/Precipitacion_animacion"
if(!dir.exists(path_prep))dir.create(path_prep)

dirs <- list.dirs("/media/miguel/TOSHIBA EXT/Biomex/")
anios <-list.files("/media/miguel/TOSHIBA EXT/Biomex/",
                   full.names = FALSE)
fils <- list.files(dirs, pattern = "bio_01.tif$", full.names = TRUE)



plan(multisession(workers = 10))


anima <- seq_along(fils)%>% furrr::future_map(function(x){
  r <- raster::raster(fils[x])
  r <- r/10
#  mexex <- extent(r)
 # mexex[2] <- -86
  pdf(paste0(path_bio1,"/",anios[x],".pdf"),
      width = 12,
      height = 10)
  pdef <- levelplot(r, main = paste0("Temperatura promedio anual ",anios[x]),
            margin=FALSE,
            colorkey=list(space="bottom"),
            par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
  print(pdef)
      dev.off()
})
plan(sequential)

plan(multisession(workers = 10))

x=1
pdefs <- list.files(path="/media/miguel/TOSHIBA EXT/Temperatura_animacion/",
                    pattern = "*.pdf",
                    full.names = TRUE)
animas<- seq_along(pdefs)%>% furrr::future_map(function(x){
    p <- pdefs[x]
    anios <- seq(1901,2016)
    rid <- image_read_pdf(p,
                          density = 300)
    pa <- "/media/miguel/TOSHIBA EXT/Animacion_temperatura_png"
    if(!dir.exists(pa))dir.create(pa)
    imapng <- magick::image_convert(image = rid,
                          format = "png")
    print(imapng)
    png(paste0(pa,"/",anios[x],".png"),
        width = 1200,
        height = 1000, res = 220)
    dev.off()
  })

  plan(sequential)
  esto <- image_join(rid)  # une imagenes
  image_animate(esto,fps=4) # se anima con una cantidad de loops
  image_write("/media/miguel/TOSHIBA EXT/Temperatura ANP´s/Figuras de temperatura por año para latex/") # se escribte en el directorio tal


################### using gifski

library(gifski)
png_files <- list.files("path/to/your/pngs/", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 1)

##################################3
