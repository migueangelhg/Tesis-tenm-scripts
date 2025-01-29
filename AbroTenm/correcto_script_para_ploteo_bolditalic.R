library(tenm)
library(raster)
library(rasterVis)
library(sf)
library(tidyselect)
library(sp)


mod_sel <- readRDS("RESULTADOS/mod_sel_tibble_abro_73occ")

layers_70_2000_dir <- ("/media/miguel/TOSHIBA EXT/Biomex_70_2000/")
m <- rgdal::readOGR("/media/miguel/TOSHIBA EXT/mapa_division_estatal_mexico/marcogeoestatal2015_gw.shp")
pabro <- rgdal::readOGR("/media/miguel/TOSHIBA EXT/abronia_IUCN_pol/data_0.shp")
a <- rio::import("abro_73occ.csv")
a <- st_as_sf(x = a,
              coords = c("decimalLongitude", "decimalLatitude"),
              crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
i <- "Abronia graminea"

mod_table <- mod_sel$mods_table
model_vars <- stringr::str_split(mod_table$fitted_vars,",")

#mod_vars <- model_vars[[47]]
mod_vars <- model_vars[[1]]
### Utilizamos el metodo predict
suit_70_2000_abro <- tenm::predict(mod_sel,model_variables=mod_vars,
                                   layers_path = layers_70_2000_dir,
                                   layers_ext = ".tif$")


######### México completo
yplot <- 5.1
xplot <- 7.2
########
#plot.new()
#par(mar = c(1,1,1,1))
pdf("RESULTADOS/abronia_proyeccion_1970_2000_mex_flat.pdf",".pdf",
    width = xplot,
    height = yplot)
########


raster::plot(suit_70_2000_abro,
             npretty = 5,
             legend.args = list(text = "Idoneidad",
                                side = 2,
                                font = 1,
                                line = .4,
                                cex = .55),
             col = rev(terrain.colors(35)), zlim = c(0, 1))

title(main = substitute(expr = bold(paste("    Distribución potencial de ",
                                          bolditalic(i)," para el período de 1970 - 2000 ")),
                        env = list(i=i)),
      cex.main = 1.05,
      font.main = 2,
      line = 1.95)

title(sub = "Idoneidad ambiental calculada a partir de las
      variables bio 04, bio 07, bio 12, bio 14 y bio 17",
      cex.sub = .48,
      font.sub = 1,
      col.sub = "black",
      line = 4.05)

title(xlab = "Longitude",
      ylab = "Latitude",
      cex.lab = .75,
      col.lab ="black",
      line = 2.6)

dev.off()

########################
##################### Plot de la región Puebla y Veracruz

yplot <- 5.97
xplot <- 4.345
pdf("RESULTADOS/abronia_proyeccion_1970_2000_iucnplot_occ_vera_pue_nosub2flat.pdf",".pdf",
    width = xplot,
    height = yplot)
raster::plot(suit_70_2000_abro,
             xlim = c(-99,-96),
             ylim = c(16, 21),
             npretty = 5,
             cex.axis=.2,
             legend.args = list(text = "Idoneidad",
                                side = 2,
                                font = 1,
                                line = .4,
                                cex = .55),
             col = rev(terrain.colors(35)), zlim = c(0, 1))
#plot(m,
#     xlim = c(-99,-96),
#     ylim = c(16, 21),
#     add = TRUE)
plot(pabro,
     border = "red",
     add = TRUE)
plot(a,
     col = "blue",
     add = TRUE,
     cex = .35,
     pch = 16)

title(main = substitute(expr = bold(paste("      Distribución potencial de ",
                                          bolditalic(i))),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.95)
title(main = substitute(expr = bold(paste(" para el período de 1970 - 2000 ")),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.15)
#title(sub = "Idoneidad ambiental calculada a partir de las
 #     variables bio 04, bio 07, bio 12, bio 14 y bio 17",
  #    cex.sub = .55,
   #   font.sub = 1,
    #  col.sub = "black",
     # line = 4.05)
title(xlab = "Longitude",
      ylab = "Latitude",
      cex.lab = .75,
      col.lab ="black",
      line = 2.6)
dev.off()









################################# Plot acercamiento a raster de idoneidad
#x11()

yplot <- 6.7
xplot <- 4.0
########
#plot.new()
#par(mar = c(0,0,0,0))
pdf("RESULTADOS/abronia_proyeccion_1970_2000_mex_pue.pdf",".pdf",
    width = xplot,
    height = yplot)
raster::plot(suit_70_2000_abro,
             xlim = c(-98,-96),
             ylim = c(17.75, 19.5),
             npretty = 5,
             legend.args = list(text = "Idoneidad",
                                side = 2,
                                font = 1,
                                line = .4,
                                cex = .55),
             col = rev(terrain.colors(35)), zlim = c(0, 1))

title(main = substitute(expr = bold(paste("      Distribución potencial de ",
                                          bolditalic(i))),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.95)

title(main = substitute(expr = bold(paste(" para el período de 1970 - 2000 ")),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.15)


#title(sub = "Idoneidad ambiental calculada a partir de las
 #     variables bio 04, bio 07, bio 12, bio 14 y bio 17",
  #    cex.sub = .55,
   #   font.sub = 1,
    #  col.sub = "black",
     # line = 4.05)

title(xlab = "Longitude",
      ylab = "Latitude",
      cex.lab = .75,
      col.lab ="black",
      line = 2.6)

dev.off()

#################### Plot de Acercamiento a poligono de experto y occs
#yplot <- 6.7
#xplot <- 4.0

yplot <- 5.5
xplot <- 3.58

pdf("RESULTADOS/abronia_proyeccion_1970_2000_iucnplot_occ_acerc_nosubflat.pdf",".pdf",
    width = xplot,
    height = yplot)
raster::plot(suit_70_2000_abro,
             cex.axis = 1.3,
             xlim = c(-97.5,-96.5),
             ylim = c(17.8, 20),
             npretty = 3,
             legend.args = list(text = "Idoneidad",
                                side = 2,
                                font = 1,
                                line = .4,
                                cex = .55),
             col = rev(terrain.colors(35)), zlim = c(0, 1))
#plot(m,
#     xlim = c(-98,-96.5),
#     ylim = c(17.8, 20),
#     add = TRUE)
plot(pabro,
     border = "red",
     add = TRUE)
plot(a,
     col = "blue",
     add = TRUE,
     cex = .45,
     pch = 16)
title(main = substitute(expr = bold(paste("      Distribución potencial de ",
                                          bolditalic(i))),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.95)
title(main = substitute(expr = bold(paste(" para el período de 1970 - 2000 ")),
                        env = list(i=i)),
      cex.main = .9,
      font.main = 2,
      line = 1.15)
#title(sub = "Idoneidad ambiental calculada a partir de las
#      variables bio 04, bio 07, bio 12, bio 14 y bio 17",
#      cex.sub = .55,
#      font.sub = 1,
#      col.sub = "black",
#      line = 4.05)
title(xlab = "Longitude",
      ylab = "Latitude",
      cex.lab = .75,
      col.lab ="black",
      line = 2.6)
dev.off()
#################################################################






















