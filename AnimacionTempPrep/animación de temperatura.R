library(dplyr)
library(raster)
install.packages("reshape2")
library(reshape2)
library(ggplot2)


##################################################
dir <- list.dirs("/media/miguel/TOSHIBA EXT/Biomex/")
fil <- list.files (dir, full.names = TRUE, pattern = "bio_01.tif$")[1:10]
seq(fil)
años <- seq(1901,1910)

r <-raster::stack(fil)
names(r) <- paste0("bio_01_",1901:1910)
### Con este codigo es posible plotear en conjunto varios rasters de varios años
crs(r)
r_df <- as.data.frame(r, xy = TRUE) %>%
  melt(id.vars = c("x","y"))

ggplot()+
  geom_raster(data = r_df, aes(x = x, y = y, fill = value))+
  facet_wrap(~variable)
#################################################################

###### Intentemos realizar lo del png de la manera
# más sencilla posible


for (x in r){
  raster::plot(xsxdvwedr3f)
}











#rp <- r %>% rasterToPoints()
rp <- data.frame(rp)


seq(rp[,-c(1,2)])
