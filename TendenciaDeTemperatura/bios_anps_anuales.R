library(raster)
library(furrr)
library(dplyr)
library(sp)

ids_pix <- rio::import("pixeles_ids_anps/pixeles_anps_ids.csv")

all_dirs <- list.dirs("E:/bioclimatic",full.names = T,recursive = F)
biomex_path <- grep(pattern = "biomex_",x = all_dirs,value = T)
x=1
plan(multisession)
bios_anps_anuales <- seq_along(biomex_path) %>% future_map_dfr(function(x){
  path_bio <- biomex_path[x]
  anio <- stringr::str_extract(string = path_bio,
                               pattern = "[0-9]+[0-9]+[0-9]+[0-9]")
  path_bio2 <- list.files(path = path_bio,
                          pattern = ".tif$",full.names = T)
  bios <- raster::stack(path_bio2)
  bios_anp <- bios[ids_pix$pixel_ID]
  df_res <- data.frame(WDPAID=ids_pix$WDPAID,anio,bios_anp)
  print(x)
  return(df_res)
  
},.progress = T)

rio::export(bios_anps_anuales,"E:/tendencia_de_temperatura_anps/bios_anps_anuales.csv")
