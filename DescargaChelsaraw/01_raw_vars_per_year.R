#-------------------------------------------------------------------------------
# Script to organize downloaded data per year
# time period: 1901-2016
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(raster)
library(stringr)
library(matrixStats)
library(furrr)
rm(list = ls())
prec_paths <- list.files("prec",pattern = ".tif$",
                         full.names = TRUE)


dirss =file.path("prec",names(n_precL))
sapply(dirss, function(x){
  if(!dir.exists(x))  dir.create(x)
})

# Data 
prec_s <- str_split(prec_paths,"_",simplify = T)

n_prec <- data.frame(year=as.numeric(prec_s[,4]),month= as.numeric(prec_s[,3]),
                     path=prec_paths)
n_precL <- n_prec %>% split(.$year)
months <- c(paste0("0",1:9),10:12)
FCOP <- seq_along(n_precL) %>% purrr::map(function(x){
  year <- names(n_precL[x])
  filess <- n_precL[[x]]$path
  file.copy(from = filess,file.path("prec",year))
  file.remove(filess)
})
                     


tmin_paths <- list.files("tmin/",pattern = ".tif$",
                         full.names = TRUE)
# Data 
tmin_s <- str_split(tmin_paths,"_",simplify = T)

n_tmin <- data.frame(year=as.numeric(tmin_s[,4]),month= as.numeric(tmin_s[,3]),
                     path=tmin_paths)

n_tminL <- n_tmin %>% split(.$year)

dirs_tmin <- file.path("tmin",names(n_tminL))
sapply(dirs_tmin, function(x){
  if(!dir.exists(x))  dir.create(x)
})



n_tminL <- n_tmin %>% split(.$year)
months <- c(paste0("0",1:9),10:12)
FCOP <- seq_along(n_tminL) %>% purrr::map(function(x){
  year <- names(n_tminL[x])
  filess <- n_tminL[[x]]$path
  file.copy(from = filess,file.path("tmin",year))
  file.remove(filess)
})





tmax_paths <- list.files("tmax/",pattern = ".tif$",
                         full.names = TRUE)
# Data 
tmax_s <- str_split(tmax_paths,"_",simplify = T)

n_tmax <- data.frame(year=as.numeric(tmax_s[,4]),month= as.numeric(tmax_s[,3]),
                     path=tmax_paths)

n_tmaxL <- n_tmax %>% split(.$year)

dirs_tmax <- file.path("tmax",names(n_tmaxL))
sapply(dirs_tmax, function(x){
  if(!dir.exists(x))  dir.create(x)
})



n_tmaxL <- n_tmax %>% split(.$year)
months <- c(paste0("0",1:9),10:12)
plan(multisession)
FCOP <- seq_along(n_tmaxL) %>% furrr::future_map(function(x){
  year <- names(n_tmaxL[x])
  filess <- n_tmaxL[[x]]$path
  file.copy(from = filess,file.path("tmax",year))
  file.remove(filess)
},.progress = T)
plan(sequential)