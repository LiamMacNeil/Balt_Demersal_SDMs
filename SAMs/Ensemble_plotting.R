library(terra)
library(raster)
library(gisland)
library(mapdata)
library(tidyverse)
library(sf)
library(oce)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Predictions
ensembles <- stack("../Predictions/Ensemble/Ensemble_BRT_RF_GAM.tif")

# Features

# Pull boundary layers (https://github.com/einarhjorleifsson/gisland)
icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  filter(subdivisio != "20") %>% 
  as("Spatial") %>% 
  crop(extent(ensembles))

'''
m <- map_data("worldHires") 

m_sf <-  st_as_sf(x = m, 
                  coords = c("long", "lat"),
                  crs = crs(icesarea)) %>% 
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

m_sf <- as(m_sf, 'Spatial') %>% 
  crop(extent(ensembles))
'''

m_sf <- read_sf("../../Oceanographic/Data/Europe_coastline_shapefile", layer = "Europe_coastline_poly")

m_sf <- m_sf %>% 
  st_transform(crs(icesarea))

m_sf <- as(m_sf, 'Spatial') %>% 
  crop(extent(icesarea))


AdultCod <- ensembles[[c("AdultCod_Q1_2001",
                         "AdultCod_Q4_2001",
                         "AdultCod_Q1_2020",
                         "AdultCod_Q4_2020")]]
JuvenileCod<- ensembles[[c(    "JuvenileCod_Q1_2001",
                               "JuvenileCod_Q4_2001",
                               "JuvenileCod_Q1_2020",
                               "JuvenileCod_Q4_2020")]]
Plaice<- ensembles[[c( "Plaice_Q1_2001",
                       "Plaice_Q4_2001",
                       "Plaice_Q1_2020",
                       "Plaice_Q4_2020")]]
Flounder<- ensembles[[c( "Flounder_Q1_2001",
                         "Flounder_Q4_2001",
                         "Flounder_Q1_2020",
                         "Flounder_Q4_2020")]]
Dab<- ensembles[[c( "Dab_Q1_2001",
                    "Dab_Q4_2001",
                    "Dab_Q1_2020",
                    "Dab_Q4_2020")]]
### Plotting

# Save collage for all years
png("../Figures/Final/Plaice_ensemble_predmap.png", width = 15, height = 12, res = 300, units = "cm")
par(mfrow=c(2,2), mar=c(0,0,3,2)+0.2)

for (i in names(Plaice)) { 
  
  rast <- Plaice[[i]]

  plot(abs(rast), 
       col=(oce.colorsViridis(20)), 
       legend.width = 1.5,
       
       #Adult cod
       #zlim = c(0.01,1.5),
       
       #Juvenile cod
       #zlim = c(0.01,1),
       
       # Dab
       #zlim = c(0.01,0.7),
       
       # Flounder
       #zlim = c(0.01,0.4),
       
       # Plaice
       zlim = c(0.01,0.4),
       
       legend.mar = 6,
       main = "",
       axes=F,
       box=F
  )
  plot(icesarea, add = T, lwd = 0.75)
  plot(m_sf, add = T, col="grey75", border = NA)
  #contour(abs((rast) > quantile(abs(rast), probs = 0.95)),  lwd = 0.175, col = "springgreen3",drawlabels = F, add=T)
  
}
dev.off()

