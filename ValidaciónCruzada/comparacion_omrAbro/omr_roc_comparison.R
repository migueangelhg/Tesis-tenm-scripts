library(tenm)
library(dplyr)
library(ntbox)

rm(list = ls())
set.seed(124)

mod_sel <- ("mod_sel_tenm_tibble_bg1.2")

lf <- list.files("/media/miguel/TOSHIBA EXT/Biomex_70_2000/", pattern = ".tif$",
                 full.names = T, recursive = F)
rs <- raster::stack(lf)

tsData <- mod_sel$temporal_df %>%  filter(year >= 1970 & year <= 2000)

cellIDs <- tenm::cells2samp(data = tsData,
                            buffer_ngbs = 100,
                            n_bg = 50000,
                            longitude = "decimalLongitude",
                            latitude = "decimalLatitude",
                            raster_mask = rs[[1]])
bg70_00<- rs[cellIDs]
bg70_00 <- data.frame(na.omit(bg))

mods_table <- mod_sel$mods_table
model_vars <- stringr::str_split(mods_table$fitted_vars,",")
vars_mod <- model_vars[[1]]


dftst <- mod_sel$temporal_df %>% filter(grepl('Test', trian_test))
dftr <- dftst <- mod_sel$temporal_df %>% filter(grepl('Train', trian_test))
cov_cent <- ntbox::cov_center(dftr[,vars_mod],
                              mve = T,
                              level = 0.975,
                              vars = vars_mod)

exstR <- raster::extract(rs,dftst[,c("decimalLongitude","decimalLatitude")])
exstR <- data.frame(exstR)
dfT70_00 <- cbind(exstR,dftst[,c("decimalLongitude","decimalLatitude", "year")])

T70_00 <- dfT70_00 %>% filter(year >= 1970 & year <= 2000)

Tn <- T70_00[,vars_mod]
#bg <- bg[,vars_mod]
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

