#-------------------------------------------------------------------------------
# Script to generate bioclimatic variables for precipitation
# time period: 1920-1980
# Uses precipitation averages from 1920 to 1980
# BIO12,BIO13,BIO14,BIO15,BIO16,BIO17
# Author: Luis Osorio-Olvera
#-------------------------------------------------------------------------------

library(furrr)
rm(list = ls())
gc()
library(magrittr)
if(!dir.exists("bioclimatic")) dir.create("bioclimatic")


prc_paths <- list.files("prec",
                       full.names = T,
                       pattern = ".tif$",
                       recursive = T)

prec_s <- str_split(prc_paths,"_",simplify = T)

n_prec <- data.frame(year=as.numeric(prec_s[,4]),month= as.numeric(prec_s[,3]),
                     path=prc_paths)
year_paths <- file.path("bioclimatic",unique(n_prec$year))
sapply(year_paths, function(x){
  if(!dir.exists(x)) dir.create(x)
})
