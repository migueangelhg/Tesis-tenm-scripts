rm(list = ls())
set.seed(444)
library(rio)
library(raster)
library(rgdal)
library(dplyr)
library(sf)
library(rgl)
library(sp)
library(ntbox)
library(ellipsenm)
library(gtools)
install.packages("Rtools4")

if (!require('devtools')) install.packages('devtools')
devtools::install_github('luismurao/ntbox')
# If you want to build vignette, install pandoc before and then
devtools::install_github('luismurao/ntbox',build_vignettes=F)

library(png)
#install.packages("ellipsenm")
#source("funciones_ayuda.R")
#Importamos los registros de presencia 
uro_cin <- rio::import("uro_occ_moransENM.csv")
#Limpieza de los duplicados espaciales por distancia
uro_cin <- ntbox::clean_dup(uro_cin,
                            longitude = "decimalLongitude",
                            latitude = "decimalLatitude",
                            threshold = .1)

###
rio::export(uro_cin, "ObjetosCsv/uroENMoccclendup.csv")
#previamente limpiados con distancia de Moran.


 

# Tranformamos nuestras presencas a objetos espaciales

sp::coordinates(uro_cin) <- ~decimalLongitude + decimalLatitude

################

am <- list.files("bioamerica_70_20", 
                 pattern = ".tif$", 
                 full.names = T,
                 recursive = T)

am <- gtools::mixedsort(am)
names(am) <- gsub("^|.tif","",am)
ames <- raster::stack(am)

#plot(ames)

# Asignamos un sistema de coordenadas de referencias
# usaremos las de las capas de worldclim
#wc1 <- raster::raster("bioamerica_70_20/bioame_bio_01.tif")
raster::crs(uro_cin) <- raster::crs(ames[[1]])
# Generamos el buffer alrededor de cada punto y unimos los poligonos

#uro_buffer <- raster::buffer(uro_occ,width =250*1000,dissolve=TRUE)
uro_buffer2 <- raster::buffer(uro_cin,width =1000*1000,dissolve=TRUE)
#uro_buffer3 <- raster::buffer(uro_occ,width =20.72964895*1000,dissolve=TRUE)

#uro_buffer_hrange <- raster::buffer(uro_cin,width =.0009772050*1000,dissolve=TRUE)
#uro_buffer_hrange <- raster::buffer(uro_cin,width =2*1000,dissolve=TRUE)
x11()
plot(ames[[1]],xlim=c(-140,-60),ylim=c(0,60))
plot(uro_buffer2,add=T,border="blue")

vars_env <- raster::mask(raster::crop(ames,uro_buffer2),uro_buffer2) %>% 
  raster::stack()

#saveRDS(vars_env,"Robjts/vars_env")
vars_env <- readRDS("Robjts/vars_env")
#vars_env <- readRDS("Robjts/vars_env")

#Colocar los puntos de entrenamiento junto con la información ambiental 
uro_cin <- as.data.frame(uro_cin)
uro_data <- raster::extract(vars_env,uro_cin[,c("decimalLongitude","decimalLatitude")])
uro_env <- data.frame(uro_cin[,c("species","decimalLongitude","decimalLatitude","year")],uro_data)
no_na_uro_env <- na.omit(data.frame(uro_env)) 

#rio::export(no_na_uro_env,"ObjetosCsv/no_na_uro_buffer_final.csv")

set.seed(444)
train_id_uro <- sample(nrow(no_na_uro_env), size = ceiling(nrow(no_na_uro_env)*0.7))

train_uro <- no_na_uro_env[train_id_uro,]
test_uro <- no_na_uro_env[-train_id_uro,]

rio::export(train_uro,"ObjetosCsv/train_uro_buffer_final.csv")
rio::export(test_uro,"ObjetosCsv/test_uro_buffer_final.csv")

##Background ambiental 

#uro_env_bg<- ntbox::sample_envbg(vars_env,nbg = 50000,rseed = 123,coordinates = TRUE)
uro_env_bg<- ntbox::sample_envbg(vars_env,nprop = 0.80,rseed = 444)
#uro_env_bg<- ntbox::sample_envbg(vars_env,nbg = 1000000,rseed = 123)

dim(uro_env_bg)
saveRDS(uro_env_bg,"Robjts/uro_buffer_bg_final")
readRDS("Robjts/uro_buffer_bg_final")
##########################
########################## Segunda parte de la modelación

no_na_uro_env <- rio::import("ObjetosCsv/no_na_uro_buffer_final.csv")
test_uro <- rio::import("ObjetosCsv/test_uro_buffer_final.csv")
train_uro <- rio::import("ObjetosCsv/train_uro_buffer_final.csv")
uro_env_bg <- readRDS("Robjts/uro_buffer_bg_final")

## Matriz de correlación 

#------------------------------------------------------------------
cor_mat <- cor(no_na_uro_env[,names(vars_env)])
env_vars <- ntbox::correlation_finder(cor_mat = cor_mat,
                                      threshold = 0.9,
                                      verbose = FALSE)$descriptors
print(env_vars)


#vars_wc <-names(uro_env_bg)#[-c(1,2)]




## Seleccion de la elipsoide 
e_selct <- ntbox::ellipsoid_selection(env_train = train_uro,
                                      env_test = test_uro,
                                      env_vars = env_vars,
                                      level = 0.975,
                                      nvarstest = c(3,4,5,6,7),
                                      env_bg = uro_env_bg,
                                      omr_criteria=0.05,
                                      parallel = TRUE,
                                      ncores = 7,
                                      comp_each = 500,
                                      proc = TRUE,sub_sample = TRUE,
                                      sub_sample_size = 500)
saveRDS(e_selct,"Robjts/e_select")
e_selct <- readRDS("Robjts/e_select")
## Filtrar los mejores modelos que pasaron las preubas de selecciÃ³n
#e_selct_omr <- e_selct %>% dplyr::filter(om_rate_train <= 0.03 & om_rate_test <= 0.03) %>% arrange(bg_prevalence)
# head(e_selct_omr)
head(e_selct)


x=8
##
#Generamos una elipse a partir de un conjunto de variables, para obtener su 
#centroide y su covarianza
bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")
vars <- unlist(bestvarcomb[[571]])
cov_centroid <- ntbox::cov_center(data = train_uro[,vars],
                                  mve = T,
                                  level = 0.975,
                                  vars = vars)
centroid <- cov_centroid$centroid
mve_cov <- cov_centroid$covariance
#options(rgl.useNULL = FALSE)
#rgl::open3d()
efit <- ntbox::ellipsoidfit2(ames[[vars]],
                             centroid = centroid,
                             covar = mve_cov,
                             plot =  TRUE,
                             level = level,
                             size = 3)
suitENMuro <- efit/raster::maxValue(efit)

x11()
plot(suitENMuro)

#options(rgl.useNULL = FALSE)
#rgl::open3d()
-
### Guardamos al modelo seleccionado 




par(mfrow=c(1,2))

rgl::plot3d(rgl::ellipse3d(cov_centroid$covariance,
                           centre = cov_centroid$centroid,
                           level = 0.99999),
            alpha=0.1,col="blue")
rgl::points3d(train_uro[,c("bioamerica_70_20.bioame_bio_01",
                           "bioamerica_70_20.bioame_bio_09",
                           "bioamerica_70_20.bioame_bio_14")])


plot3d(efit)







#############################################################################
bestvarcomb <- stringr::str_split(e_selct$fitted_vars,",")
modelos <- length(bestvarcomb)%>% purrr::map(function(x){
  vars1 <- unlist(bestvarcomb[[x]])
  
  cov_centroid <- ntbox::cov_center(train_uro[,vars1],
                                    mve = T,level = 0.975,
                                    vars = vars1 )
  centroid <- cov_centroid$centroid
  
  mve_cov <- cov_centroid$covariance
  options(rgl.useNULL = FALSE)
  rgl::open3d()
  efit <- ntbox::ellipsoidfit2(ames[[vars1]],
                               centroid = centroid,
                               covar = mve_cov,
                               plot =  TRUE,
                               level = level,size = 3)
  plot(efit)
  return(efit)
  
})

x11()
plot(modelos)

raster::writeRaster(efit, "E:/urocyon_niche/Final_Uro_wc_70_2000/raster_idoneidad.tif")
#rio::export(efit,"E:/urocyon_niche/Final_Uro_wc_70_2000/idoneidadx_1_bufferfinal.csv")
modelos <- raster::stack(modelos)
windows()
raster::plot(modelos$Suitability.1)


raster::plot(modelos$Suitability)
plot()
windows()
raster::plot(efit$Suitability)