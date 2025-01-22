setwd("E:/bioclimatic")

getwd()
mex <- shapefile("E:/mex_gadm/gadm36_MEX_0.shp")
dirs_names <- list.dirs(".",full.names = F,recursive = F)
list.dirs(".",full.names = F,recursive=F)
dirs_names <- list.dirs(".",full.names = F,recursive=F)
recortes <- seq_along(dirs_names) %>% purrr::map(function(x){
  dirnombre <- paste0("biomex_",dirs_names[x])
  dir.create(dirnombre)
  capas_path <- list.files(dirs_names[x],full.names = T,pattern = "tif$")
  capas_stack <- raster::stack(capas_path)
  biomex_stack <- raster::mask(raster::crop(capas_stack,mex),mex)
  print("recorte hecho")
  nombres_raster <- file.path(dirnombre,paste0("biomex_",names(capas_stack),".tif"))
  raster::writeRaster(biomex_stack,filename = nombres_raster, bylayer =T, overwrite=T)
  return()
})
x=1
?writeRaster
