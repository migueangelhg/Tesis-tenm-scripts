library(hsi)
library(dplyr)
library(rgl)
rm(list = ls())
source("scripts/00_functions.R")
dd <- rio::import("occ_Urocyon_cinereoargenteus.csv")
datos_sucios <- nrow(dd) 
d3 <- dd %>% filter(year >= 1970 & year <= 2000)
datos_en_perido_wc2 <- nrow(d3)

d4 <- ntbox::clean_dup(data = d3,longitude = "decimalLongitude", 
                       latitude =  "decimalLatitude",
                       threshold = (0.5/60)*2)

datos_curados <- nrow(d4)

nrow(d4)
d2 <- dd %>% select(decimalLongitude,decimalLatitude,year)
test_sp <- hsi::sp_temporal_data(occs=d2,longitude = "decimalLongitude",
                                 latitude = "decimalLatitude",
                                 sp_year_var="year",
                                 layers_by_year_dir ="bioclimatic",
                                 layers_ext = "*.tif$",
                                 reclass_year_data = F)
test_sp_clean <- hsi::clean_dup_by_year(this_species = test_sp,
                                        threshold = (0.5/60)*2)

nrow(test_sp_clean$oocs_data)
#d3 <- test_sp_clean$oocs_data
#hsi::sp_temporal_data(occs=d3,longitude = "decimalLongitude",
#                      latitude = "decimalLatitude",
#                      sp_year_var="year",
#                      layers_by_year_dir ="bioclimatic",
#                      layers_ext = "*.tif$",
#                      reclass_year_data = F)
sp_logfile <- "Urocyon_cinereoargenteus.log"
sp_name <- "Urocyon_cinereoargenteus"
cat("El numero de registros sin duplicados espacio-temporales para",sp_name,
    "es",nrow(test_sp_clean$sp_coords),"\n",file = sp_logfile,append = T)
cat("Distribucion de presencias por anio:",file = sp_logfile,append = T)
capture.output(table(test_sp_clean$sp_occs_year),file = sp_logfile,append = T)
e_test <- ex_by_year(this_species = test_sp_clean,
                     layers_pattern = ".tif$")

env_bg <- bg_by_year(this_species = e_test,
                     layers_pattern = ".tif$",
                     ncores = 10,
                     buffer_ngbs = 10,
                     n_bg = 60000)


occ_env <- na.omit(e_test$coords_env_data_all[,-c(1:3,ncol(e_test$coords_env_data_all))])
mat_cor <- cor(occ_env)
vars_cor <- ntbox::correlation_finder(mat_cor,threshold = 0.8,verbose = F)
vars_cor$descriptors
vars_cor$list_cor

seleccion_elp <- ntbox::ellipsoid_selection(env_train = e_test$env_data_train,
                                            env_test = e_test$env_data_test,
                                            env_vars = vars_cor$descriptors,
                                            nvarstest = c(3),
                                            level = 0.975,
                                            mve = T,
                                            env_bg = env_bg,
                                            omr_criteria = 0.05,
                                            parallel = TRUE,
                                            ncores = 20,
                                            comp_each = 300,
                                            proc = TRUE,proc_iter = 200)

head(seleccion_elp)

vars_mods <- stringr::str_split(seleccion_elp$fitted_vars,pattern = ",")
mod_nicho <- cov_center(na.omit(e_test$env_data_train[,c("bio_01",vars_mods[[4]])]),
                        mve = T,level = 0.975,1:3)


plot3d(env_bg[,vars_mods[[4]]])
points3d(e_test$coords_env_data_all[,vars_mods[[4]]],col="red")
rgl::plot3d(rgl::ellipse3d(mod_nicho$covariance,centre=mod_nicho$centroid,
                           level = 0.999),
            add=T,alpha=0.3,col="blue")
points3d(e_test$coords_env_data_all[,vars_mods[[4]]],col="red")
points3d(e_test$coords_env_data_all[,vars_mods[[4]]],col="red")

bios_path <- list.files("bioclimatic/2016",pattern = ".tif$",full.names = T)
bios <- raster::stack(bios_path)
bios <- bios[[vars_mods[[4]]]]

modelo_nicho <- ntbox::ellipsoidfit2(bios,centroid = mod_nicho$centroid,
                                     covar = mod_nicho$covariance,level = 0.9999,
                                     plot = TRUE,size = 3)

raster::writeRaster(modelo_nicho,filename = "modelo_de_juguete.tif",overwrite=T)

e_test$env_data_test
tiempo_seleccion <- system.time({
  best_model_2004 <- find_best_model_ntbox2(this_species = e_test,
                                                cor_threshold = 0.8,
                                                nbg_points = 50000,
                                                omr_criteria = 0.2,
                                                ellipsoid_level = 0.975,
                                                nvars_to_fit = 3,
                                                E = 5,
                                                RandomPercent = 70,
                                                NoOfIteration = 50,
                                                parallel=T,
                                                n_cores=20)
})
