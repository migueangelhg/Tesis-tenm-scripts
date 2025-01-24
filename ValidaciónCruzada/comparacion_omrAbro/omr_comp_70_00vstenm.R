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
AbrDataCl <- rio::import("obj_csv/Abronia_occ_gbif_doi_clean_dup.csv")

cellIDs <- tenm::cells2samp(data = AbrDataCl,
                            buffer_ngbs = 100,
                            n_bg = 50000,
                            longitude = "decimalLongitude",
                            latitude = "decimalLatitude",
                            raster_mask = rs[[1]])
bg70_00<- rs[cellIDs]
bg70_00 <- data.frame(na.omit(bg70_00))


bestvarcomb <- stringr::str_split(e_sel$fitted_vars,",")
vars_mod <- unlist(bestvarcomb[[2]])

dftst <- rio::import("obj_csv/test_abro.csv")
dftr <- rio::import("obj_csv/train_abro.csv")

cov_cent <- ntbox::cov_center(dftr[,vars_mod],
                              mve = T,
                              level = 0.975,
                              vars = vars_mod)

exstR <- raster::extract(rs,dftst[,c("decimalLongitude","decimalLatitude")])
exstR <- data.frame(exstR)
dfT70_00 <- cbind(exstR,dftst[,c("decimalLongitude","decimalLatitude", "year")])

T70_00 <- dfT70_00 %>% filter(year >= 1970 & year <= 2000)

Tn <- T70_00[,vars_mod]

bg70_00 <- bg70_00[,vars_mod]
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



abro1901_2016occ <- rio::import("abro_modsel_clean.csv")

ex <- raster::extract(rs,abro1901_2016occ[,c("decimalLongitude","decimalLatitude")])
exdf <- data.frame(ex)
dfex <- cbind(exdf,abro1901_2016occ[,c("decimalLongitude","decimalLatitude", "year")])

bgtenm <- dfex[,vars_mod]

suit_bgtenm <- exp(-0.5 * mahalanobis(bgtenm,
                                  center = cov_cent$centroid,
                                  cov = cov_cent$covariance))

suit_tst70_00 <- exp(-0.5 * mahalanobis(Tn,
                                   center = cov_cent$centroid,
                                   cov = cov_cent$covariance))


rocp <- ntbox::pROC(continuous_mod = suit_bgtenm,
                    test_data = suit_tst70_00,
                    n_iter = 1000,
                    E_percent = 5,
                    boost_percent = 50,
                    rseed = T,
                    sub_sample = T,
                    sub_sample_size = 10000)

rocp$pROC_summary

E_slt$in_Ellipsoid
