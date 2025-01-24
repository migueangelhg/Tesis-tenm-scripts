rm(list = ls())

library(tenm)
library(dplyr)
library(ntbox)
# Establesco una semilla para alinear el reloj interno
set.seed(124)
#Primeramente estoy asignando al objeto e_select  con el archivo donde
#se encuentran los modelos candidatos con AUC y Omr, del cual la tengo
#previamente seleccionado el modelo correcto ENM.
e_sel <- readRDS("e_select")

## Es es el directorio a las capas con las que se construyo el modelo
#de 1970 a 2000. (ENM)
lf <- list.files("D:/bioamerica_70_20/", pattern = ".tif$",
                 full.names = T, recursive = F)
lf <- gtools::mixedsort(lf)
ames <- raster::stack(lf)
lf <- gsub("^|.tif","",lf)
lf <- gsub("D:/", "", lf)
lf <- gsub("//","/", lf)
#lf <- gsub("biamerica_70_20/")

#nomb <- paste0("bioamerica_70_20.",gsub(".tif", "", lf))
#Se realizo un conjunto (stack) raster con las 19 capas para poder manipularlas
#ames <- raster::stack(lf)
names(ames) <- lf
###################3
### Todos estos datos de ocurrencia ya estaban previamente filtrados del 70
# al 2000
#Base de datos de las presencias totales con el clean dup de ENM
Ur70_20DataCl <- rio::import("D:/Miguelon/Urocyon/tradicional/UroENM70_2000/ObjetosCsv/uroENMoccclendup.csv")
#Buffer a partir todas las presencias limpias de U.cin
cellIDs <- tenm::cells2samp(data = Ur70_20DataCl,
                            buffer_ngbs = 100,
                            n_bg = 200000,
                            longitude = "decimalLongitude",
                            latitude = "decimalLatitude",
                            raster_mask = ames[[1]])

#ID pixel de las capas de america
bg70_00 <- ames[cellIDs]
#Eliminamos todas aquellas celdas sin información ambiental
bg70_00 <- data.frame(na.omit(bg70_00))
#Elegimos las variables para el mejor modelo de ENM

bestvarcomb <- stringr::str_split(e_sel$fitted_vars,",")
vars_modENM <- unlist(bestvarcomb[[571]])
#Objetos de test del ENM
dftst <- rio::import("D:/Miguelon/Urocyon/tradicional/UroENM70_2000/ObjetosCsv/test_uro_buffer_final.csv")
#oBJETOS DE TRAIN
dftr <- rio::import("D://Miguelon/Urocyon/tradicional/UroENM70_2000/ObjetosCsv/train_uro_buffer_final.csv")

cov_cent <- ntbox::cov_center(dftr[,vars_modENM],
                              mve = T,
                              level = 0.975,
                              vars = vars_modENM)
#Extraemos la informacion ambiental de los datos de TEST
exstR <- raster::extract(ames,dftst[,c("decimalLongitude","decimalLatitude")])
exstR <- data.frame(exstR)
dfT70_00 <- cbind(exstR,dftst[,c("decimalLongitude","decimalLatitude", "year")])

#Info ambiental con años y georeferencias
T70_00 <- dfT70_00 %>% filter(year >= 1970 & year <= 2000)

#Unicamente información ambiental de los registros Train
Tn <- T70_00[,vars_modENM]

bg70_00 <- bg70_00[,vars_modENM]
suit_bg <- exp(-0.5 * mahalanobis(bg70_00,
                                  center = cov_cent$centroid,
                                  cov = cov_cent$covariance))

suit_tst <- exp(-0.5 * mahalanobis(Tn,
                                   center = cov_cent$centroid,
                                   cov = cov_cent$covariance))


rocp <- ntbox::pROC(continuous_mod = suit_bg,
                    test_data = suit_tst,
                    n_iter = 1000,
                    E_percent = 5,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 10000)

rocp$pROC_summary

E_slt <- ntbox::inEllipsoid(centroid = cov_cent$centroid,
                            eShape = cov_cent$covariance,
                            env_data = Tn,
                            level = 0.975)

E_slt$in_Ellipsoid


######TENM
#Objeto R data con todos los modelos seleccionados para TENM
mod_sel <- readRDS("mod_sel_uro")
# Directorio a las capas chelsa 1901 a 2016

layersTS_dir <- ("D:/bioamerica/")
#Obtencin de los registros de presencia
uro1901_2016occClean <- mod_sel$temporal_df

uroTenmBg <- mod_sel$env_bg



#dfOcctenm <- data.frame(uro1901_2016occClean[,c("decimalLongitude",
 #                                    "decimalLatitude",
  #                                   "year")])
#tenmOcc<-tenm::sp_temporal_data(dfOcctenm,
 #                      longitude = "decimalLongitude",
  #                     latitude = "decimalLatitude",
   #                    sp_date_var = "year",
    #                   occ_date_format = "y",
     #                  layers_date_format = "y",
      #                 layers_by_date_dir = layersTS_dir,
       #                layers_ext = "*.tif")

#tenm::bg_by_date(tenmOcc,
 #                buffer_ngbs = 10,
  #               n_bg =1000)
#Filtrado de los registros tanto de test como de train
urtraintenm <- uro1901_2016occClean %>% filter(grepl('Train', trian_test))
urtsttenm <- uro1901_2016occClean %>% filter(grepl('Test', trian_test))
#Elegimos el mejor modelo (#1)
mods_table <- mod_sel$mods_table
bestvarcomb <- stringr::str_split(mods_table$fitted_vars,",")
# Vector de caracter con los nombres de las variables
# del mejor modelo tenm
vars_modTENM <-unlist(bestvarcomb[[1]])
#bg 1901_2016
exTrnTS <- data.frame(urtraintenm[,vars_modTENM])
exTstTS <- data.frame(urtsttenm[,vars_modTENM])
BGtenm <- uroTenmBg[,vars_modTENM]

# Calculamos el centroide del nicho a partir de las
# variables para el mejor modelo de TENM
cov_centTenm <- ntbox::cov_center(na.omit(urtraintenm[,vars_modTENM]),
                              mve = T,
                              level = 0.975,
                              vars = vars_modTENM)



#extstTenm <- raster::extract(ames,urtsttenm[,c("decimalLongitude","decimalLatitude")])
#extstdf <- data.frame(na.omit(extstTenm))
#dfextst <- cbind(extstTenm,urtsttenm[,c("decimalLongitude","decimalLatitude", "year")])

#dfextst <- data.frame(na.omit(dfextst))

#dfextst <- dfextst %>% filter(year >= 1970 & year <= 2000)

#vars_modTENM <- paste0("bioamerica_70_20.", vars_modTENM)

#bgtenm <- dfextst[,vars_modTENM]

suit_bgtenm <- exp(-0.5 * mahalanobis(BGtenm,
                                      center = cov_centTenm$centroid,
                                      cov = cov_centTenm$covariance))

suit_tstTenm <- exp(-0.5 * mahalanobis(exTstTS,
                                            center = cov_centTenm$centroid,
                                            cov = cov_centTenm$covariance))


rocptenm <- ntbox::pROC(continuous_mod = suit_bgtenm,
                    test_data = suit_tstTenm,
                    n_iter = 100,
                    E_percent = 5,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 100)

rocptenm$pROC_summary

E_sltTenm <- ntbox::inEllipsoid(centroid = cov_centTenm$centroid,
                            eShape = cov_centTenm$covariance,
                            env_data = na.omit(exTstTS),
                            level = 0.975)

E_sltTenm$in_Ellipsoid
E_sltTenm <- ntbox::inEllipsoid(centroid = cov_centTenm$centroid,
                                eShape = cov_centTenm$covariance,
                                env_data = dftst[,c(5,6,12)],
                                level = 0.975)
