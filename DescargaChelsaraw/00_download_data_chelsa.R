#-------------------------------------------------------------------------------
# Script to download time series data of CHELSA
# Precipitation, tmin and tmax 
# time period: 1901-2016
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(raster)
library(stringr)
library(furrr)
rm(list = ls())
gc()
prec_info <- readLines("envidatS3paths.txt")
prec_info <- gsub(" ","",prec_info)

prec_years <- grep("prec",
                   prec_info,value = TRUE)

prec_split <- str_split(prec_years,pattern = "/",simplify = T)
prec_names <- prec_split [,ncol(prec_split)]

if(!dir.exists("prec")) dir.create("prec")
prec_paths <- file.path("prec",prec_names)

data_prec <- data.frame(url =prec_years, varname =prec_names,varpath =prec_paths)

plan(multisession)
x=50
prec_down <- 1:nrow(data_prec) %>% furrr::future_map(function(x){
  download.file(data_prec$url[x],
                destfile = data_prec$varpath[x],
                method = "wget")
  return()
},.progress = TRUE)
plan(sequential)



tmax_info <- readLines("envidatS3paths.txt")
tmax_info <- gsub(" ","",tmax_info)

tmax_years <- grep("tmax",
                   tmax_info,value = TRUE)

tmax_split <- str_split(tmax_years,pattern = "/",simplify = T)
tmax_names <- tmax_split [,ncol(tmax_split)]

if(!dir.exists("tmax")) dir.create("tmax")
tmax_paths <- file.path("tmax",tmax_names)
data_tmax <- data.frame(url =tmax_years, 
                        varname =tmax_names,
                        varpath =tmax_paths)
plan(multisession)
tmax_down <- 1:nrow(data_tmax) %>% furrr::future_map(function(x){
  download.file(data_tmax$url[x],
                destfile = data_tmax$varpath[x],
                method = "wget")
  return()
},.progress = TRUE)
plan(setequal)


tmin_info <- readLines("envidatS3paths.txt")
tmin_info <- gsub(" ","",tmin_info)

tmin_years <- grep("tmin",
                   tmin_info,value = TRUE)

tmin_split <- str_split(tmin_years,pattern = "/",simplify = T)
tmin_names <- tmin_split [,ncol(tmin_split)]

if(!dir.exists("tmin")) dir.create("tmin")
tmin_paths <- file.path("tmin",tmin_names)
data_tmin <- data.frame(url =tmin_years, 
                        varname =tmin_names,
                        varpath =tmin_paths)
plan(multisession)
tmin_down <- 1:nrow(data_tmin) %>% furrr::future_map(function(x){
  download.file(data_tmin$url[x],
                destfile = data_tmin$varpath[x],
                method = "wget")
  return()
},.progress = TRUE)
plan(setequal)
