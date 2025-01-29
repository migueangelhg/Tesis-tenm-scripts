library(furrr)
#?dismo::biovars
rm(list = ls())
gc()
library(magrittr)

prc_path <- list.files("prec_avg_1920_1980",
                       full.names = T,pattern = ".tif$")
tmn_path <- list.files("tmin_avg_1920_1980",
                       full.names = T,pattern = ".tif$")

rtmp <- raster::raster(tmn_path[1])
rprc <- raster::raster(prc_path[1])

rtmp_vals <- rtmp[]
rprc_vals <- rprc[]
nona_tmp <- which(!is.na(rtmp_vals))
nona_prc <- which(!is.na(rprc_vals))
nonas <- intersect(nona_tmp,nona_prc)
#nona_mat <- as.matrix(nonas)
#rio::export(nona_mat,"nona_cells.rds")
rm(c("rtmp","rprc",
     "rtmp_vals","rprc_vals",
     "nona_tmp","nona_prc"))
gc()

plan(multisession(workers = 2))
options(future.globals.maxSize = 8500 * 1024^2)
prc <- seq_along(prc_path) %>% furrr::future_map_dfc(function(x){
  pre <- raster::raster(prc_path[x])
  pre_vals <- pre[]
  pre_vals <- pre_vals[nonas]
  df1 <- data.frame(pre_vals)
  names(df1) <- paste0("pre",x)
  return(df1)
},.progress = TRUE)
plan(sequential)
prc <- as.matrix(prc)

x <- 1:12
xmax = max(x)
cuartos <- lapply(seq_along(x), function(i){
  if(i==1) return(x[1:3])
  if(i==(xmax-1)) return(c(x[i],x[i+1],x[1]))
  if(i==xmax) return(c(x[i],x[1],x[2]))
  else return(c(x[i],x[i+1],x[i+2]))
})


prc_por_cuarto <- seq_along(cuartos) %>% purrr::map_dfc(function(x){
  r00 <- matrixStats::rowSums2(prc[,cuartos[[x]]])
  prc_cuarto <- data.frame(r00)
  names(prc_cuarto) <- paste0(cuartos[[x]],collapse = "_")
  print(x)
  return(prc_cuarto)
})
rm(list = c("prc"))
gc()
prc_por_cuarto <- as.matrix(prc_por_cuarto)


ids_max <- Rfast::rowMaxs(prc_por_cuarto)
ids_min <- Rfast::rowMins(prc_por_cuarto)
prc_cuartos_id <- data.frame(ids_max,ids_min,nonas)
prc_cuartos_id <- as.matrix(prc_cuartos_id)
mode(prc_cuartos_id) <- "integer"
rio::export(prc_cuartos_id,"prc_cuartos_ids_1920_1980.rds")



