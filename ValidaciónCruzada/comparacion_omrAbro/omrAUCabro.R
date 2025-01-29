rm(list = ls())

library(tenm)
library(dplyr)
library(ntbox)

set.seed(124)

e_sel <- readRDS("e_slect_abro")

lf <- list.files("D:/Biomex_70_2000/", pattern = ".tif$",
                 full.names = T, recursive = F)
rs <- raster::stack(lf)

### Todos estos datos de ocurrencia ya estaban previamente filtrados del 70
# al 2000
AbrDataClENM <- rio::import("obj_csv/Abronia_occ_gbif_doi_clean_dup.csv")

#cellIDs <- tenm::cells2samp(data = AbrDataClENM,
  #                          buffer_ngbs = 100,
   #                         n_bg = 1000,
    #                        longitude = "decimalLongitude",
     #                       latitude = "decimalLatitude",
      #                      raster_mask = rs[[1]])
#bg70_00<- rs[cellIDs]
#bg70_00 <- data.frame(na.omit(bg70_00))

bg70_00 <- rio::import("obj_csv/back_abro_10km.csv")

bestvarcomb <- stringr::str_split(e_sel$fitted_vars,",")
vars_modENM <- unlist(bestvarcomb[[2]])

dftst <- rio::import("obj_csv/test_abro.csv")
dftr <- rio::import("obj_csv/train_abro.csv")

cov_cent <- ntbox::cov_center(dftr[,vars_modENM],
                              mve = T,
                              level = 0.975,
                              vars = vars_modENM)

exstR <- raster::extract(rs,dftst[,c("decimalLongitude","decimalLatitude")])
exstR <- data.frame(exstR)
dfT70_00 <- cbind(exstR,dftst[,c("decimalLongitude","decimalLatitude", "year")])



Tn <- dfT70_00[,vars_modENM]

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
                    E_percent = 1,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 1000)

rocp$pROC_summary

E_sltENM <- ntbox::inEllipsoid(centroid = cov_cent$centroid,
                            eShape = cov_cent$covariance,
                            env_data = Tn,
                            level = 0.975)

E_sltENM$in_Ellipsoid


### Lista de puntos de presencia Abro
#Llamamos al objeto mod_sel

mod_sel <- readRDS("mod_sel_tenm_tibble_bg1.2")

mod_table <- mod_sel$mods_table
bestvarcombTENM <- stringr::str_split(mod_table$fitted_vars,",")
vars_modTENM <- unlist(bestvarcombTENM[[203]])

dftenm <- data.frame(mod_sel$temporal_df)
BGtenm <- data.frame(mod_sel$env_bg)

dftrnts <- data.frame(dftenm %>% filter(grepl('Train', trian_test)))
dftstts <- data.frame(dftenm %>% filter(grepl('Test', trian_test)))
#abro1901_2016occ <- rio::import("abro_modsel_clean.csv")
#ex <- raster::extract(rs,abro1901_2016occ[,c("decimalLongitude","decimalLatitude")])
#exdf <- data.frame(ex)
#dfex <- cbind(exdf,abro1901_2016occ[,c("decimalLongitude","decimalLatitude", "year")])

#bgtenm <- dfex[,vars_mod]

bgtenm <- BGtenm[,vars_modTENM]
TstenvTenm <- dftstts[,vars_modTENM]
## Calibración del modelo con los puntos de train
cov_centTENM <- ntbox::cov_center(dftrnts[,vars_modTENM],
                              mve = T,
                              level = 0.975,
                              vars = vars_modTENM)



suit_bgtenm <- exp(-0.5 * mahalanobis(bgtenm,
                                      center = cov_centTENM$centroid,
                                      cov = cov_centTENM$covariance))

suit_tsttenm <- exp(-0.5 * mahalanobis(TstenvTenm,
                                        center = cov_centTENM$centroid,
                                        cov = cov_centTENM$covariance))

#Calculo de la curva ROC a partir de la prevalencia
rocpTenm <- ntbox::pROC(continuous_mod = suit_bgtenm,
                    test_data = suit_tsttenm,
                    n_iter = 10000,
                    E_percent = 1,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 100000)

rocpTenm$pROC_summary

E_sltTENM <- ntbox::inEllipsoid(centroid = cov_centTENM$centroid,
                               eShape = cov_centTENM$covariance,
                               env_data = TstenvTenm,
                               level = 0.975)

E_sltTENM$in_Ellipsoid

##################### Comparativo #############################
#Centroide tenm con puntos test de ENM#
E_sltENMvsTenm <- ntbox::inEllipsoid(centroid = cov_centTENM$centroid,
                                     eShape = cov_centTENM$covariance,
                                     env_data = dftst[,vars_modTENM],
                                     level = 0.975)

omrEnmvstenm <- E_sltENMvsTenm$in_Ellipsoid
length(which(omrEnmvstenm %in% 0))/length(omrEnmvstenm)
#centroide ENM con puntos test de TENM de varsENM
#exTstTS <- data.frame(urtsttenm[,c("bioame_bio_03","bioame_bio_09","bioame_bio_15")])
E_sltTenmvsENM <- ntbox::inEllipsoid(centroid = cov_cent$centroid,
                                     eShape = cov_cent$covariance,
                                     env_data = dftstts[,vars_modENM],
                                     level = 0.975)
om <-E_sltTenmvsENM$in_Ellipsoid

dfFin <- cbind(dftstts,E_sltTenmvsENM$in_Ellipsoid)

length(which(om %in% 0))/length(om)

