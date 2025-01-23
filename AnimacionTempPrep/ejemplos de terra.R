#Ejemplos del paquete terra terra::plot

f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
plot(r)
e <- c(6.37, 6.41, 49.9, 50.1)
plot(r, plg=list(ext=e, title="Title\n", title.cex=1.25),
     pax=list(side=1:4, retro=TRUE))
north("topleft")

##Raster examples


r <- raster(nrows=10, ncol=10)
r <- setValues(r, 1:ncell(r))
r <- raster::raster(r)
plot(r)

############

r <- raster(nrows=10, ncols=10)
r <- setValues(r, 1:ncell(r))
plot(r)

e <- extent(r)
plot(e, add=TRUE, col='red', lwd=4)
e <- e / 2
plot(e, add=TRUE, col='red')

r2 <- sqrt(r)
plot(r, r2)
plot(r, r2, gridded=TRUE)
