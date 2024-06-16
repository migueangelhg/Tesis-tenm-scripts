library(raster)
library(dplyr)
install.packages("rasterVis")
library(rasterVis)
install.packages("Cairo")
library(Cairo)
library(ggplot2)

path_bio1 <- "/media/miguel/TOSHIBA EXT/Temperatura_animacion"
if(!dir.exists(path_bio1))dir.create(path_bio1)
path_prep <- "/media/miguel/TOSHIBA EXT/Precipitacion_animacion"
if(!dir.exists(path_prep))dir.create(path_prep)

dirs <- list.dirs("/media/miguel/TOSHIBA EXT/Biomex/")
anios <-list.files("/media/miguel/TOSHIBA EXT/Biomex/",
                   full.names = FALSE)
fils <- list.files(dirs, pattern = "bio_01.tif$", full.names = TRUE)
#r <- raster::stack(fils)
#años <- seq(1901:2016)

#names(r) <-paste0("bio_01_",1901:1903)
x=1
anima <- seq_along(fils)%>% furrr::future_map(function(x){
  r <- raster::raster(fils[x])
  r <- r/10
  mexex <- extent(r)
  mexex[2] <- -86
  reso <- 1200
  length <- 3.25*reso/72
  #png(".png",units="in",res=reso,height=length,width=length)
  #par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
  #plot(r,avg="vertical",spread.estimate="stddev",col="black",lty=3, lwd=3)
  #dev.off()
  #length <- 1200
  par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
  raster::plot(r,
               legend.args = list(text = "Temperatura (°C)",
                                  side = 4,
                                  font = 2,
                                  line = 2.5,
                                  cex = 0.8))
  png(paste0(path_bio1,"/",anios[x],".png"),
      width = length,
      height = length, res = 220)
  dev.off()

  levelplot(r, main = paste0("Temperatura promedio anual ",anios[x]),
            margin=FALSE,
            colorkey=list(space="bottom"),
            par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
  #plot(mexex, col = ("transparent"), xlab = "", ylab  = "")
  #plot(r, xlim = c(mexex[1],-82),
   #    add =TRUE,
    #   main = paste0("Temperatura promedio anual ",anios[x]),
     #  horizontal = TRUE,
      # legend.args = list(text = "°C", side = 3,
       #                   font = 2, line = .3, cex = 0.8))


  dev.off()

})


####
ggsave(
  "/media/miguel/TOSHIBA EXT/Temperatura_animacion/",
  r,
  width = 3.25,
  height = 3.25,
  dpi = 1200
###

#raster::plot(r, legend.args = list(text = "Temperatura (°C)", side = 4,
               #    font = 2, line = 2.5, cex = 0.8))
#r %>% raster::plot(r[1],main = paste0("Promedio de temperatura anual"),años)
#png("",width = 1200, height = 1000, res = 220,)

#raster::plot(r[1],
 #            xmn =
  #           x = "longitude",
   #          y = "latitude")
