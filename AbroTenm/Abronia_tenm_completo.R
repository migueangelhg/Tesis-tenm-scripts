#Para  comenzar llamemos las librerias que vamos a utilizar
library(tenm)
library(tidyverse)
library(rio)
library(dplyr)
library(raster)

#Primeramente importamos los registros de A.graminea
#Después de haber sido curados.
#Agregamos una semilla para alinear el reloj interno
set.seed(124)
#Archivo csv
#abro <- rio::import("Abronia occ/abronia_graminea_sin quintana roo.csv")
#abro <- abro %>% filter(year >= 1901, year <= 2016)
# rio::export("Abronia occ/abro_gram_1901_2016.csv")
#abro <- rio::import("Abronia occ/abro_gram_1901_2016.csv")
abro <- rio::import("Abronia occ/abro_gram_1901_2016_noelevmax_nonas.csv")

d <- geodata::elevation_3s(lon = abro$decimalLongitude[1],
                           lat = abro$decimalLatitude[1],
                           path = "/media/miguel/TOSHIBA EXT/A_graminea_tenm/Abronia_tenm_3,4,5,6,7 modsel/elevabro.tif")
el_abro_elev <- raster::extract(d,abro[c("decimalLongitude","decimalLatitude")])
abroelevdfadd <- cbind(abro, el_abro_elev)
abro_fil <- dplyr::filter(abroelevdfadd, srtm_17_09 >= 1200 & srtm_17_09 <= 2900)
abro_fil <- abro_fil %>% dplyr::filter(!is.na("srtm_17_09"))
# Eliminamos los NA
abro <-  abro_fil %>% filter(!is.na("decimalLatitude"))
#Colocamos opcionalmente el directorio completo de las capas en un objeto
laers_dir <- "/media/miguel/TOSHIBA EXT/Biomex/"
# Convertimos esos registros de presencia en un objeto sp_temporal_data
#Indicamos el nombre de las columnas que contienen la información de long y lat
abrostd <- tenm::sp_temporal_data(occs = abro,
                                  longitude = "decimalLongitude",
                                  latitude = "decimalLatitude",
                                  sp_date_var = "year",
                                  occ_date_format = "y",
                                  layers_date_format = "y",
                                  layers_by_date_dir = laers_dir,
                                  layers_ext = "*tif$")

## Posteriormente limpiamos los duplicados espaciales por fecha.
#Agregamos el raster .tif de la resolucion a la que queremos limpiar
#abrostd_clean <- tenm::clean_dup_by_date(abrostd, threshold = 0.55555/60)
#10/60)
#Agregamos el raster .tif de la resolucion a la que queremos limpiar
rmasc <- raster::raster("/media/miguel/TOSHIBA EXT/Biomex/1901/biomex_bio_01.tif")
#Limpiamos duplicados, 0 significa en el mismo pixel
abrostd_clean <- tenm::clean_dup_by_date(this_species = abrostd,
                                          by_mask = TRUE,
                                          n_ngbs = 0,
                                          raster_mask = rmasc)
saveRDS(abrostd_clean, file = "/media/miguel/TOSHIBA EXT/A_graminea_tenm/bgExtract70_2000/abrostd_clean")

abrodf <- as.data.frame(abrostd_clean$temporal_df)
#0.55555/60 significa 0.009259167 grados, es decir 1.029619 de diametro y por
#tanto 1.617322 km2 cuadrados de area

## Vamos a proceder a extraer la información ambiental
# con un proceso en parelelo.
future::plan("multisession", workers = 12)
abroex <- tenm::ex_by_date(this_species = abrostd_clean,
                           train_prop = 0.7)
future::plan("sequential")
#Correlación de variables
varscors <- tenm::correlation_finder(environmental_data =
                                     abroex$env_data[,-ncol(abroex$env_data)],
                                     method = "spearman",
                                     threshold = 0.8,
                                     verbose = FALSE)
vars2fit <- varscors$descriptors
#vars2fit3 <- c("biomex_bio_02","biomex_bio_04","biomex_bio_03","biomex_bio_07",
 #              "biomex_bio_12","biomex_bio_13","biomex_bio_14","biomex_bio_15",
  #             "biomex_bio_16","biomex_bio_17","biomex_bio_18","biomex_bio_19")
### Obtenemos el background ambiental en paralelo
future::plan("multisession", workers = 12)
abbg <- tenm::bg_by_date(abroex,
                         buffer_ngbs = 36,
                         n_bg = 100000)
future::plan("sequential")

#Selección de modelos.
mod_sel <- tenm::tenm_selection(this_species = abbg,
                                omr_criteria = 0.1,
                                ellipsoid_level = 0.975,
                                vars2fit = vars2fit,
                                nvars_to_fit =c(3,4,5,6,7),
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration = 1000,
                                parallel = TRUE,
                                n_cores = 12)
saveRDS(mod_sel,
        file = "RESULTADOS/objetosR_tenmfinal/mod_sel_tenm_tibble")


#Simplemente guardamos en un objeto el directorio de los archivos de la información ambiental
layers_70_2000_dir <- ("/media/miguel/TOSHIBA EXT/Biomex_70_2000/")
### Utilizamos el metodo predict
suit_70_2000 <- tenm::predict(mod_sel,model_variables=NULL,
                              layers_path = layers_70_2000_dir,
                              layers_ext = ".tif$")
plot(suit_70_2000)

### Guardar el archivo RASTER de ideneoidad-suitability

raster::writeRaster(suit_70_2000, "RESULTADOS/Abronia_idoneidad_70_2000_buffer_mod_1000_3ngbs_73occ.tif", format = "tif")

###Archivo shapefile
mx_es <- rgdal::readOGR("/media/miguel/TOSHIBA EXT/mapa_division_estatal_mexico/marcogeoestatal2015_gw.shp")

#### Acercamiento
m_es <- rgdal::readOGR("/media/miguel/TOSHIBA EXT/mapa_division_estatal_mexico/marcogeoestatal2015_gw.shp")
xplot2 <- 4.5
yplot2 <- 5.1

#x11()

plot.new()
pdf("abronia_proyeccion_1970_2000_acercamiento_mod.pdf",".pdf",
    width = xplot2,
    height = yplot2)

rst <- raster::plot(suit_70_2000,
                    main = "Distribución potencial de Abronia graminea
para el periodo de 1970-2000",
                    cex.main = 0.9,
                    cex.lab = .8,
                    cex.axis = .5,
                    xlim = c(-99, -95),
                    ylim = c(16, 21),
                    xlab = "Longitude",
                    ylab = "Latitude",
                    legend.args = list(text = "Idoneidad",
                                       side = 2,
                                       font = 1,
                                       line = .4,
                                       cex = .8),
                    maxpixels=500000,
                    npretty = 5,
                    cex = 10,
                  #  col.regions = col,
                   # colorkey = list(width = .6, height = .6),
                    col = rev(terrain.colors(35)), zlim = c(0, 1))#,
#add = TRUE)

mx <- sp::plot(m_es,
           xlim = c(-99, -95),
           ylim = c(16, 21),
           cex = 4,
           add = TRUE)
print(mx)#,
           #add = TRUE,)
dev.off()














###plot de la infromacion
xplot <- 7.7
yplot <- 4.91
x11()
plot.new()
pdf("abronia_proyeccion_1970_2000_threshold_mod.pdf",".pdf",
    width = xplot,
    height = yplot)

rst <- raster::plot(suit_70_2000,
                     main = title("Distribución potencial de Abronia graminea para el periodo de 1970-2000"),
                     xlim = c(-108, -85),
                     ylim = c(14, 25),
                     xlab = "Longitude",
                     ylab = "Latitude",
                     legend.args = list(text = "Idoneidad",
                                        side = 2,
                                        font = 1,
                                        line = .45,
                                        cex = .8),
                     maxpixels=50000000,
                     npretty = 10,
                     cex = 2,
                     col = rev(terrain.colors(35)))#,
#add = TRUE)
print(rst)
mx <- plot(mx_es, add = TRUE)
dev.off()
