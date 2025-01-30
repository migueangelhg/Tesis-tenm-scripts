library(tenm)
library(future)
library(rio)
library(tibble)
library(raster)
rm(list = ls())
set.seed(111)
#importamos la infromacion de los registros de presencia.
uro <- rio::import("uro_occ_tenm_MornsDecade.csv")
#um <- rio::import("Limpieza estructurada por morans y años/uro_occ_tenm_modelf.csv")
# Directorio donde se encuentran las carpetas por año para las 19 bioclim
dir_layers <-  "/media/miguel/TOSHIBA EXT/bioamerica"
### Raster de la mascara para la limpieza
rmasc <-raster::raster("/media/miguel//TOSHIBA EXT/bioamerica/1901/bioame_bio_01.tif")

#Transformamos los registros en objetos de la clase sp.temporal.modeling.
urostd <- tenm::sp_temporal_data(occs = uro,
                                 longitude = "decimalLongitude",
                                 latitude = "decimalLatitude",
                                 sp_date_var = "year",
                                 occ_date_format = "y",
                                 layers_date_format = "y",
                                 layers_by_date_dir = dir_layers,
                                 layers_ext = "*.tif$")

#Limpieza de duplicados espaciales por fecha

urost_clean <- tenm::clean_dup_by_date(this_species = urostd,
                                       n_ngbs = 50,
                                       by_mask = TRUE,
                                       raster_mask = rmasc)

# Extracción de la información ambiental
future::plan("multisession", workers = 12)
ubex <- tenm::ex_by_date(this_species = urost_clean,
                         train_prop = 0.7)
future::plan("sequential")

saveRDS(ubex,
        file = "ObjetosRData/uro_extract")
ubex <- readRDS("ObjetosRData/uro_extract.gz")
## Correlacion de variables.
varcors <- tenm::correlation_finder(environmental_data = ubex$env_data[,-ncol(ubex$env_data)],
                                    method = "spearman",
                                    threshold = 0.9,
                                    verbose = FALSE)
vars2fit <- varcors$descriptors
#background ambiental siempre siendo cuidadosos de no exceder el
#numero de nucleos o workers para que la memoria RAM no se vea rebasada.

future::plan("multisession", workers = 12)
ubbg <- tenm::bg_by_date(ubex,buffer_ngbs = 800, n_bg = 10000)
future::plan("sequential")


#
#rio::export(ubbg,"/media/miguel/TOSHIBA EXT/Urocyon_tenm_final/(tenm) Urocyon_noviembre/", ".csv")

### La selección de modelos
mod_sel <- tenm::tenm_selection(this_species = ubbg,
                                omr_criteria = 0.1,
                                ellipsoid_level = 0.975,
                                vars2fit = vars2fit,
                                nvars_to_fit =c(3,4,5,6,7),
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration = 1000,
                                parallel = TRUE,
                                n_cores = 16)


dat_abro <- as.data.frame(mod_sel)


write.table(mod_sel$temporal_df, file = "tem_df_tibble_mod_sel.csv")
write.table(mod_sel$mods_table, file = "mod_sel_df.csv")
write.table(mod_sel$env_bg, file = "env_bg_df.csv")


rio::export("urocion_mod_selec.csv",
            format = ".csv")
layers_70_2000_dir <- ("/media/ecoinflab/TOSHIBA EXT/bioamerica_70_20/")

suit_70_2000 <- tenm::predict(mod_sel,model_variables=NULL,
                              layers_path = layers_70_2000_dir,
                              layers_ext = ".tif$")


mod_table <- mod_sel$mods_table
model_vars <- stringr::str_split(mod_table$fitted_vars,",")


suit_70_2000_mod_sel_160 <- tenm::predict(mod_sel,model_variables=model_vars[[24]],
                                          layers_path = layers_70_2000_dir,
                                          layers_ext = ".tif$")


rgl::rglwidget()
x11()
plot(suit_70_2000_mod_sel_160)
#raster::writeRaster(suit_70_2000_mod_sel_160, "urocyon_tenm_70_00_rank2_14_09_05_bio.tif")



plot(suit_70_2000)
rgl::rglwidget()
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_rank1_14_08_05_bio.tif")





suit_70_2000_mod_sel_10 <- tenm::predict(mod_sel,model_variables=model_vars[[16]],
                                         layers_path = layers_70_2000_dir,
                                         layers_ext = ".tif$")
rgl::rglwidget()
plot(suit_70_2000_mod_sel_10)

raster::writeRastbackground ambiental siempre siendo cuidadosos de no exceder el
#numero de nucleos o workers para que la memoria RAM no se vea rebasada.

future::plan("multisession", workers = 12)
ubbg <- tenm::bg_by_date(ubex,buffer_ngbs = 20, n_bg = 50000)
future::plan("sequential")


#
#rio::export(ubbg,"/media/miguel/TOSHIBA EXT/Urocyon_tenm_final/(tenm) Urocyon_noviembre/", ".csv")

### La selección de modelos
mod_sel <- tenm::tenm_selection(this_species = ubbg,
                                omr_criteria = 0.1,
                                ellipsoid_level = 0.975,
                                vars2fit = vars2fit,
                                nvars_to_fit =c(3,4,5),
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration = 1000,
                                parallel = TRUE,
                                n_cores = 16)


dat_abro <- as.data.frame(mod_sel)


write.table(mod_sel$temporal_df, file = "tem_df_tibble_mod_sel.csv")
write.table(mod_sel$mods_table, file = "mod_sel_df.csv")
write.table(mod_sel$env_bg, file = "env_bg_df.csv")


rio::export("urocion_mod_selec.csv",
            format = ".csv")
layers_70_2000_dir <- ("/media/ecoinflab/TOSHIBA EXT/bioamerica_70_20/")

suit_70_2000 <- tenm::predict(mod_sel,model_variables=NULL,
                              layers_path = layers_70_2000_dir,
                              layers_ext = ".tif$")


mod_table <- mod_sel$mods_table
model_vars <- stringr::str_split(mod_table$fitted_vars,",")


suit_70_2000_mod_sel_160 <- tenm::predict(mod_sel,model_variables=model_vars[[2]],
                                          layers_path = layers_70_2000_dir,
                                          layers_ext = ".tif$")


rgl::rglwidget()

plot(suit_70_2000_mod_sel_160)
#raster::writeRaster(suit_70_2000_mod_sel_160, "urocyon_tenm_70_00_rank2_14_09_05_bio.tif")



plot(suit_70_2000)
rgl::rglwidget()
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_rank1_14_08_05_bio.tif")





suit_70_2000_mod_sel_10 <- tenm::predict(mod_sel,model_variables=model_vars[[16]],
                                         layers_path = layers_70_2000_dir,
                                         layers_ext = ".tif$")
rgl::rglwidget()
plot(suit_70_2000_mod_sel_10)

raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_rank10_14_07__09_01_bio.tif")






urorastenm <- raster::plot(suit_70_2000, xlim = c(-175, -25), ylim = c(0,60))
print(urorastenm)
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_pred_nat.tif", ".tif")
.#############
library(rasterVis)
png("urocyon_proyeccion70_00(5000bg)_levelplot.png",width = 1200, height = 1000, res = 220)
pdef <- levelplot(suit_70_2000, main = "Distribución potencial de U.cinereoargenteus
                  para el periodo 1970-2000",
                  margin = FALSE,
                  colorkey=list(space="bottom"),
                  par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
print(pdef)
dev.off()
############

png("urocyon_proyeccion70_00(5000bg)2.png",width = 1200, height = 1000, res = 220)
raster::plot(suit_70_2000)
dev.off()
er(suit_70_2000, "urocyon_tenm_70_00_rank10_14_07__09_01_bio.tif")






urorastenm <- raster::plot(suit_70_2000, xlim = c(-175, -25), ylim = c(0,60))
print(urorastenm)
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_pred_nat.tif", ".tif")
.#############
library(rasterVis)
png("urocyon_proyeccion70_00(5000bg)_levelplot.png",width = 1200, height = 1000, res = 220)
pdef <- levelplot(suit_70_2000, main = "Distribución potencial de U.cinereoargenteus
                  para el periodo 1970-2000",
                  margin = FALSE,
                  colorkey=list(space="bottom"),
                  par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
print(pdef)
dev.off()
############

png("urocyon_proyeccion70_00(5000bg)2.png",width = 1200, height = 1000, res = 220)
raster::plot(suit_70_2000)
dev.off()
rs$descriptors

### Obtenemos el background ambiental siempre siendo cuidadosos de no exceder el
#numero de nucleos o workers para que la memoria RAM no se vea rebasada.

future::plan("multisession", workers = 12)
ubbg <- tenm::bg_by_date(ubex,buffer_ngbs = 20, n_bg = 50000)
future::plan("sequential")


#
#rio::export(ubbg,"/media/miguel/TOSHIBA EXT/Urocyon_tenm_final/(tenm) Urocyon_noviembre/", ".csv")

### La selección de modelos
mod_sel <- tenm::tenm_selection(this_species = ubbg,
                                omr_criteria = 0.1,
                                ellipsoid_level = 0.975,
                                vars2fit = vars2fit,
                                nvars_to_fit =c(3,4,5),
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration = 1000,
                                parallel = TRUE,
                                n_cores = 16)


dat_abro <- as.data.frame(mod_sel)


write.table(mod_sel$temporal_df, file = "tem_df_tibble_mod_sel.csv")
write.table(mod_sel$mods_table, file = "mod_sel_df.csv")
write.table(mod_sel$env_bg, file = "env_bg_df.csv")


rio::export("urocion_mod_selec.csv",
            format = ".csv")
layers_70_2000_dir <- ("/media/ecoinflab/TOSHIBA EXT/bioamerica_70_20/")

suit_70_2000 <- tenm::predict(mod_sel,model_variables=NULL,
                              layers_path = layers_70_2000_dir,
                              layers_ext = ".tif$")


mod_table <- mod_sel$mods_table
model_vars <- stringr::str_split(mod_table$fitted_vars,",")


suit_70_2000_mod_sel_160 <- tenm::predict(mod_sel,model_variables=model_vars[[2]],
                              layers_path = layers_70_2000_dir,
                              layers_ext = ".tif$")


rgl::rglwidget()

plot(suit_70_2000_mod_sel_160)
#raster::writeRaster(suit_70_2000_mod_sel_160, "urocyon_tenm_70_00_rank2_14_09_05_bio.tif")



plot(suit_70_2000)
rgl::rglwidget()
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_rank1_14_08_05_bio.tif")





suit_70_2000_mod_sel_10 <- tenm::predict(mod_sel,model_variables=model_vars[[16]],
                                          layers_path = layers_70_2000_dir,
                                          layers_ext = ".tif$")
rgl::rglwidget()
plot(suit_70_2000_mod_sel_10)

raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_rank10_14_07__09_01_bio.tif")






urorastenm <- raster::plot(suit_70_2000, xlim = c(-175, -25), ylim = c(0,60))
print(urorastenm)
raster::writeRaster(suit_70_2000, "urocyon_tenm_70_00_pred_nat.tif", ".tif")
.#############
library(rasterVis)
png("urocyon_proyeccion70_00(5000bg)_levelplot.png",width = 1200, height = 1000, res = 220)
pdef <- levelplot(suit_70_2000, main = "Distribución potencial de U.cinereoargenteus
                  para el periodo 1970-2000",
                  margin = FALSE,
                  colorkey=list(space="bottom"),
                  par.settings = RdBuTheme(region = rev(RColorBrewer::brewer.pal(9, 'RdBu'))))
print(pdef)
dev.off()
############

png("urocyon_proyeccion70_00(5000bg)2.png",width = 1200, height = 1000, res = 220)
raster::plot(suit_70_2000)
dev.off()
