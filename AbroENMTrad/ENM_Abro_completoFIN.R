rm(list = ls())

library(tenm)
library(ntbox)
library(raster)
library(dplyr)

set.seed(666)

#Informacion ambiental 
mx <- list.files("D:/Biomex_70_2000/",
                 pattern = ".tif$",
                 full.names = TRUE,
                 recursive = FALSE)
mx_ch <- raster::stack(mx)

#Genero un objeto que contenga los registros de presencia
abro_occ <- rio::import("Abro_occ/Abronia_occ_gbif_doi.csv")

# Procedemos a generar otro objeto de abronia graminea, 
#filtrados de  1970 a 2000 
abro70_00 <- abro_occ %>% dplyr::filter(year >=1970 & year <= 2000)

#ELiminamos todas las filas con presencia de NA en las coordenadas de latitud
abro <-  abro70_00 %>% dplyr::filter(!is.na(abro70_00$decimalLatitude))


# Realizamos un filtrado deacuerdo a la elevacion 
# Primeramente necesitamos obtener el raster de elevacion con respecto 
# a una presencia y guardalo en algun directorio 

#d <- geodata::elevation_3s(lon = abro$decimalLongitude[1],
 #                          lat = abro$decimalLatitude[1],
  #                         path = "Elevation_raster/elevras.tif/")

### Si ya está descargado podemos solo llamarlo 
d <- raster::raster("Elevation_raster/elevras.tif/srtm_17_09.tif")

srtm_17_09 <- raster::extract(d,abro[c("decimalLongitude","decimalLatitude")])
#Unimos la informaci?n de  elevaci?n con el df original 
abro_occ <- cbind(abro, srtm_17_09)

#Filtramos las coordenas con elevaciones entre el rango de 1200 y 2900 msnm
abro_fil <- dplyr::filter(abro_occ, srtm_17_09 >= 1200 & srtm_17_09 <= 2900)
#quitamos a las NA de elevacion
abro_occ <- abro_fil %>% dplyr::filter(!is.na(abro_fil$srtm_17_09))

abro_cl <- tenm::clean_dup(data = abro_occ,
                longitude = "decimalLongitude",
                latitude = "decimalLatitude",
                by_mask = TRUE,
                raster_mask = raster::raster(mx[[1]]),
                n_ngbs = 0)

#Volvemos los registros en objetos espaciales 
# Coordenadas 
sp::coordinates(abro_cl) <- ~decimalLongitude + decimalLatitude

########### Buffer ambiental ###############

# Asignamos un sistema de coordenadas de referencias
# usaremos las de las capas promedio como sistema de coordenadas
raster::crs(abro_cl) <- raster::crs(mx_ch[[1]])
 
# Generamos el buffer alrededor de cada punto y unimos los poligonos

abro_buffer <- raster::buffer(abro_cl, width = 60*1000,dissolve=TRUE)
dim(abro_buffer)
# Observamos el plot del buffer
#windows()
x11()
plot(mx_ch[[1]],xlim=c(-105,-90),ylim=c(15,25))
plot(abro_buffer,add=T,border="blue")
#Variables ambientales 
vars_env <- raster::mask(raster::crop(mx_ch,abro_buffer),abro_buffer) %>%
  raster::stack()

#Colocar los registros de presencia con la informacion ambiental 
#Convertir el abro_occc a dataframe en lugar de SpatialPointsDataFrame
abro_occ <- as.data.frame(abro_cl)
abro_data <- raster::extract(mx_ch,abro_occ[,c("decimalLongitude","decimalLatitude")])
abro_env <- data.frame(abro_occ[,c("decimalLongitude","decimalLatitude","year")],abro_data)
no_na_abro_env <- na.omit(data.frame(abro_env))

rio::export(no_na_abro_env,"objetos/abro_env_no_na.csv")


#Colocar los puntos de entrenamiento 
#set.seed(123)
train_id_abro <- sample(nrow(no_na_abro_env), size = ceiling(nrow(no_na_abro_env)*0.7))
train_abro <- no_na_abro_env[train_id_abro,]
test_abro <- no_na_abro_env[-train_id_abro,]
rio::export(train_abro,"objetos/train_abro.csv")
rio::export(test_abro,"objetos/test_abro.csv")
#Background ambiental 
#abro_env_bg <- ntbox::sample_envbg(vars_env,nbg = 50000, rseed = 666,coordinates = TRUE)
abro_env_bg <- ntbox::sample_envbg(vars_env,nprop = .9, rseed = 666,coordinates = TRUE)

dim(abro_env_bg)
#rio::export(abro_env_bg,"objetos/back_abro_120km.csv")

########################## Segunda parte de la modelacion ##########################

no_na_abro_env <- rio::import("objetos/abro_env_no_na.csv")
test_abro <- rio::import("objetos/test_abro.csv")
train_abro <- rio::import("objetos/train_abro.csv")
abro_env_bg <- rio::import("objetos/back_abro_120km.csv")

cor_mat <- cor(no_na_abro_env[,names(mx_ch)])

env_vars <- ntbox::correlation_finder(cor_mat = cor_mat,
                                      threshold = 0.9,
                                      verbose = FALSE)$descriptors
print(env_vars)


## Selección de elipsoide 

e_selct <- ntbox::ellipsoid_selection(env_train = train_abro,
                                      env_test = test_abro,
                                      env_vars = env_vars,
                                      level = 0.975,
                                      nvarstest = c(2,3,4),
                                      env_bg = abro_env_bg,
                                      omr_criteria=0.1,
                                      parallel = F,
                                      ncores = 9,
                                      comp_each = 500,
                                      proc = TRUE, 
                                      proc_iter = 500,
                                      sub_sample = TRUE,
                                      sub_sample_size = 1000)
##############################################3
View(e_selct)
########################################################################


## Filtrar los mejores modelos que pasaron las preubas de selecciÃ³n
#e_selct_omr <- e_selct %>% dplyr::filter(om_rate_train <= 0.2 & om_rate_test <= 0.2) %>% arrange(bg_prevalence)

saveRDS(e_selct, "e_slect_abro")
e <-readRDS("e_slect_abro")
View(e)

bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")
vars1 <- unlist(bestvarcomb[[2]])
cov_centroid <- ntbox::cov_center(train_abro[vars1],
                                  mve = T,level = 0.975,
                                  vars = vars1)
centroid <- cov_centroid$centroid
mve_cov <- cov_centroid$covariance
#options(rgl.useNull = FALSE)
#rgl::open3d()
efit <- ntbox::ellipsoidfit2(mx_ch[[vars1]],
                             centroid = centroid,
                             covar = mve_cov,
                             plot = TRUE,
                             level = level, size = 3)
x11()
plot(efit,
     xlim = c(-97.5,-96.75),
     ylim = c(17.75, 19.5))



x11()
plot(efit,
     xlim = c(-99.5,-96),
     ylim = c(17, 20))
