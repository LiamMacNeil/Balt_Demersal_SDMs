library(raster)
library(rgdal)
library(prioritizr)
library(prioritizrdata)
library(scales)
library(stars)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  as("Spatial")

Cod <- raster("../BO_Cod/Static/Outputs/Predicted_Static.tif")

Velocity <- raster("VoCC/voccMag.tif")
Velocity <- crop(mask(Velocity, icesarea), extent(icesarea))
Velocitytresampled <- projectRaster(Velocity,Cod,method = 'bilinear')

grid <- rasterToPolygons(Velocitytresampled)
  
Balt <- read_sf("../Balt_shape" , layer = "Baltic")
pr <- as(Balt, "Spatial")
plot(pr)
MPAs <- read_sf("../Balt_MPAs", layer = "HELCOM_MPAs_2019_2")

MPAs <- st_transform(MPAs, crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")

MPArast <-as(st_rasterize(MPAs), "Raster")

# Save MPA raster layers containing continuous values
for (p in 1:nlayers(MPArast)) {
  writeRaster(MPArast[[p]], filename = paste0(names(MPArast)[p], ".tif"))
}

MPArast_resamp <- crop(mask(MPArast, icesarea), extent(Velocitytresampled))
