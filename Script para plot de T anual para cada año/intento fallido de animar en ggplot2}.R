#install.packages("ecmwfr")
library("ecmwfr")
library("raster")
library("tidyverse")
#install.packages("rnaturalearth")
library("rnaturalearth")
library("sf")
#install.packages("gganimate")
library("gganimate")
library("ggplot2")
#install.packages("gifski")
library("gifski")

###########################################

#library(showtext)
#font_add_google("Prata", regular.wt = 400)
#showtext_auto()

#request <- list(
#  date = "2021-07-21/2022-07-28",
#  type = "forecast",
#  format = "netcdf_zip",
#  variable = "black_carbon_aerosol_optical_depth_550nm",
#  time = "00:00",
#  leadtime_hour = as.character(1:120),
#  area = c(90, -180, 0, 180),
#  dataset_short_name = "cams-global-atmospheric-composition-forecasts",
#  target = "download.zip"
#)


#file <- wf_request(request,
#                   user = "316282666@iztacala.unam.mx")


###########################################3

añoschr <-as.character(seq(1901,2016))
añosdirs <- list.dirs("/media/miguel/TOSHIBA EXT/Biomex/")

añospaths <- list.files(añosdirs, full.names = TRUE, pattern = "bio_01.tif$")
bios <- raster::stack(añospaths)
#dirs <- lapply(años, FUN = list.files)
robinson <- CRS("+proj=robin +over")
#Create a bounding box for the robinson proyection

bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin=-160,
            xmax=-70,
            ymax=45,
            ymin=10),crs = st_crs(4326)),n=100))

bb_robinson <- st_transform(bb, as.character(robinson))


tem <- raster::stack(bios)

tem_df <- tem %>%
  projectRaster(., res = 50000, crs = robinson)%>%
  rasterToPoints%>%
  as.data.frame()



mex_map <- geom_tile(
    data = tem_df,
    aes(
      x = x,
      y = y,
      fill = log(val),
      group = date)) +
  scale_fill_viridis_c(
    option = "B"
  ) +
  geom_sf(data=bb_robinson,
        colour='white',
        linetype='solid',
        fill = NA,
        size= 1)
