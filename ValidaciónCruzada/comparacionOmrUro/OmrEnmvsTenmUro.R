#OmrEnmvsTenmUro#
rm(list = ls())
library(tenm)
library(dplyr)
library(ntbox)
set.seed(124)

#Objeto con los modelos finales de ENM tradiconales
e_sel <- readRDS("e_select")
lf <- list.files("D:/bioamerica_70_20/", pattern = ".tif$",
                 full.names = T, recursive = F)
rs <- raster::stack(lf)
lf <- list.files("D:/bioamerica_70_20/", pattern = ".tif$",
                 full.names = T, recursive = F)
lf <- gtools::mixedsort(lf)
ames <- raster::stack(lf)
lf <- gsub("^|.tif","",lf)
lf <- gsub("D:/", "", lf)
lf <- gsub("//","/", lf)
names(rs) <- lf
### Todos estos datos de ocurrencia ya estaban previamente filtrados del 70
# al 2000

#Puntos de presencia ENM urocyon
UroDataClENM <- rio::import("D:/Miguelon/Urocyon/tradicional/UroENM70_2000/ObjetosCsv/uroENMoccclendup.csv")
#cellIDs <- tenm::cells2samp(data = AbrDataClENM,
#                          buffer_ngbs = 100,
#                         n_bg = 1000,
#                        longitude = "decimalLongitude",
#                       latitude = "decimalLatitude",
#                      raster_mask = rs[[1]])
#bg70_00<- rs[cellIDs]
#bg70_00 <- data.frame(na.omit(bg70_00))

#Objeto con lel BG de los modelos finales.
bg70_00ENMUro <- readRDS("D:/Miguelon/Urocyon/tradicional/UroENM70_2000/Robjts/uro_buffer_bg_final")
bestvarcomb <- stringr::str_split(e_sel$fitted_vars,",")

#Combinacion de variables del modelo final para los ENM

vars_modENM <- unlist(bestvarcomb[[571]])

#de 1970 a 2000. (ENM)
#Particion de datos train (calibracion) y test (validacion)
dftstUroENM <- rio::import("ObjtCsv/test_uro_buffer_final.csv")
dftrUroENM <- rio::import("ObjtCsv/train_uro_buffer_final.csv")

#Calculo de la elipsoide de volumen minima con datos train
# ENM tradicionales
cov_cent <- ntbox::cov_center(dftrUroENM[,vars_modENM],
                              mve = T,
                              level = 0.975,
                              vars = vars_modENM)
# Extraccion de informaion ambiental
exstRuro <- raster::extract(rs,dftstUroENM[,c("decimalLongitude","decimalLatitude")])
exstRuro <- data.frame(exstRuro)
# datos test
dfT70_00enm <- cbind(exstRuro,dftstUroENM[,c("decimalLongitude","decimalLatitude", "year")])
#Informacion ambiental datos test
Tn <- dfT70_00enm[,vars_modENM]

#Informacion ambiental del BG
bg70_00ENMUro <- bg70_00ENMUro[,vars_modENM]

#Calculo de idoneidad BG ENM
suit_bg <- exp(-0.5 * mahalanobis(bg70_00ENMUro,
                                  center = cov_cent$centroid,
                                  cov = cov_cent$covariance))

#Calculo de idoneidad test ENM
suit_tst <- exp(-0.5 * mahalanobis(Tn,
                                   center = cov_cent$centroid,
                                   cov = cov_cent$covariance))
# Calculo de la curva ROC para ENM utilizando en el parametro de
# "continuous_mod" la idoneidad del BG anteriormente importado para
# ENM, mientras que como "test_data" al modelo de idoneidad con
# datos de validacion test ENM

rocp <- ntbox::pROC(continuous_mod = suit_bg,
                    test_data = suit_tst,
                    n_iter = 1000,
                    E_percent = 1,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 1000)
rocp$pROC_summary

#Calculo de la distancia de los puntos Test ENM tradicional
# al centroide del modelo final
E_sltENM <- ntbox::inEllipsoid(centroid = cov_cent$centroid,
                               eShape = cov_cent$covariance,
                               env_data = Tn,
                               level = 0.975)

#Cantidad de datos dentro de la elipse de volumen minimo

omrENM <- E_sltENM$in_Ellipsoid

# Calculo de la omisión utilizando la proporcion de
#datos fuera de la elipse entre el total de todos los puntos

length(which(omrENM %in% 0))/length(omrENM)

##################### procedimiento tenm ####################
#Llamamos al objeto mod_sel (tenm)
mod_sel <- readRDS("mod_sel_uro")
mod_table <- mod_sel$mods_table
bestvarcombTENM <- stringr::str_split(mod_table$fitted_vars,",")
#Combinación de variables del modelo final tenm

vars_modTENM <- unlist(bestvarcombTENM[[1]])

#Extraccion de las presencias o registros tenm a patir del objeto
#"mod_sel_uro" que las contiene

dftenm <- data.frame(mod_sel$temporal_df)

#extracción del BG tenm

BGtenm <- data.frame(mod_sel$env_bg)

#Partición de los datos de calibracion (train) y validacion (test) en tenm

dftrnts <- data.frame(dftenm %>% filter(grepl('Train', trian_test)))
dftstts <- data.frame(na.omit(dftenm %>% filter(grepl('Test', trian_test))))
#abro1901_2016occ <- rio::import("abro_modsel_clean.csv")
#ex <- raster::extract(rs,abro1901_2016occ[,c("decimalLongitude","decimalLatitude")])
#exdf <- data.frame(ex)
#dfex <- cbind(exdf,abro1901_2016occ[,c("decimalLongitude","decimalLatitude", "year")])
#bgtenm <- dfex[,vars_mod]

#Bg de la combinacion de variables del modelo final

bgtenm <- BGtenm[,vars_modTENM]

#informacion ambiental de los puntos test

TstenvTenm <- dftstts[,vars_modTENM]

## Calibración del centroide del modelo tenm con los puntos de train
cov_centTENM <- ntbox::cov_center(na.omit(dftrnts[,vars_modTENM]),
                                  mve = T,
                                  level = 0.975,
                                  vars = vars_modTENM)
#Calculo de idoneidad BG tenm
suit_bgtenm <- exp(-0.5 * mahalanobis(bgtenm,
                                      center = cov_centTENM$centroid,
                                      cov = cov_centTENM$covariance))
#Calculo de idoneidad test tenm
suit_tsttenm <- exp(-0.5 * mahalanobis(TstenvTenm,
                                       center = cov_centTENM$centroid,
                                       cov = cov_centTENM$covariance))
#calculo de la curva roc para el algoritmo tenm
#Utilizando como modelo continuo (continuous_mod) el BG
#y como test_data el modeo continuo con puntos de validación tenm
rocpTenm <- ntbox::pROC(continuous_mod = suit_bgtenm,
                        test_data = suit_tsttenm,
                        n_iter = 10000,
                        E_percent = 1,
                        boost_percent = 50,
                        rseed = T,
                        sub_sample = T,
                        sub_sample_size = 100000)
rocpTenm$pROC_summary
#Calculo de las distancias al centro de la elipse (TENM)
#para determinar si se encuentra dentro o fuera de ella
E_sltTENM <- ntbox::inEllipsoid(centroid = cov_centTENM$centroid,
                                eShape = cov_centTENM$covariance,
                                env_data = na.omit(TstenvTenm),
                                level = 0.975)
E_sltTENM$in_Ellipsoid
omrTENM <- E_sltENM$in_Ellipsoid
#Omision: cantidad de datos de validacion (test) presentes dentro
# elipsede volumen minimo calculda con una proporcion de 0,975 de
# los datos de Train (calibracion)
length(which(E_sltTENM %in% 0))/length(omrTENM)


##################### Comparativo #############################
#Centroide tenm con puntos test de ENM#
varsTenmtoVs <- paste0("bioamerica_70_20.",vars_modTENM)
E_sltENMvsTenm <- ntbox::inEllipsoid(centroid = cov_centTENM$centroid,
                                     eShape = cov_centTENM$covariance,
                                     env_data = dftstUroENM[,varsTenmtoVs],
                                     level = 0.975)

omrTENMvsenm <- E_sltENMvsTenm$in_Ellipsoid
length(which(omrTENMvsenm %in% 0))/length(omrTENMvsenm)
#centroide ENM con puntos test de TENM de varsENM
#exTstTS <- data.frame(urtsttenm[,c("bioame_bio_03","bioame_bio_09","bioame_bio_15")])
varsENMvstenm <- gsub("bioamerica_70_20.", "",vars_modENM)
E_sltTenmvsENM <- ntbox::inEllipsoid(centroid = cov_cent$centroid,
                                     eShape = cov_cent$covariance,
                                     env_data = na.omit(dftstts[,varsENMvstenm]),
                                     level = 0.975)
om <-E_sltTenmvsENM$in_Ellipsoid
dfFin <- cbind(dftstts,E_sltTenmvsENM$in_Ellipsoid)
length(which(om %in% 0))/length(om)

