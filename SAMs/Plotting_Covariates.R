library(terra)
library(raster)
library(oce)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 2001-2010
Q1_layers_cut <- stack(list.files("../Predictors/Monthly/2001-2010/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))
Q4_layers_cut <- stack(list.files("../Predictors/Monthly/2001-2010/Q4/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))

#2011-2020
Q1_layers_cut <- stack(list.files("../Predictors/Monthly/2011-2020/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))
Q4_layers_cut <- stack(list.files("../Predictors/Monthly/2011-2020/Q4/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))


names(Q1_layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                          "MLD" ,"Seabed.Habitat","Surface.oxygen" , "Surface.salinity","Surface.temperature")

names(Q4_layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                          "MLD" ,"Seabed.Habitat","Surface.oxygen" , "Surface.salinity","Surface.temperature")

# Extra (crude) plotting script to reproduce 

###
# Bottom Oxygen
###
png("../Predictors/BotOxygen_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[1]], col = oce.colorsOxygen(20), zlim = c(-4,13),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()
png("../Predictors/BotOxygen_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[1]], col = oce.colorsOxygen(20), zlim = c(-4,13),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Surface Oxygen
###
png("../Predictors/SurOxygen_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[7]], col = oce.colorsOxygen(20), zlim = c(5,15),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()
png("../Predictors/SurOxygen_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[7]], col = oce.colorsOxygen(20), zlim = c(5,15),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Surface Temp
###
png("../Predictors/SurTemp_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[9]], col = oce.colorsTemperature(20), zlim = c(0,12),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

png("../Predictors/SurTemp_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[9]], col = oce.colorsTemperature(20), zlim = c(0,12),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Bottom Temp
###
png("../Predictors/BotTemp_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[3]], col = oce.colorsTemperature(20), zlim = c(0,14),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

png("../Predictors/BotTemp_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[3]], col = oce.colorsTemperature(20), zlim = c(0,14),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Surface Salin
###
png("../Predictors/SurSalin_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[8]], col = oce.colorsSalinity(20), zlim = c(5,34.5),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

png("../Predictors/SurSalin_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[8]], col = oce.colorsSalinity(20), zlim = c(5,34.5),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Bottom Salin
###
png("../Predictors/BotSalin_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[2]], col = oce.colorsSalinity(20), zlim = c(5,36),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

png("../Predictors/BotSalin_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[2]], col = oce.colorsSalinity(20), zlim = c(5,36),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

###
# Chloropyhll
###
png("../Predictors/Chl_Q1_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q1_layers_cut[[4]], col = oce.colorsChlorophyll(20), zlim = c(0,4),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()

png("../Predictors/Chl_Q4_2020.png", width = 16, height = 15, units = "cm", res = 300)
plot(Q4_layers_cut[[4]], col = oce.colorsChlorophyll(20), zlim = c(0,4),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75),colNA = "white")
dev.off()






##########
# Extra salinity plot for entire time period
##########
sob <- stack(list.files("../Predictors/Monthly/2011-2020/sob", 
                        pattern = "\\.tif$", full.names = TRUE))
sob <- mean(sob)



r <- stack(list.files("../Predictors/Monthly/2001-2010/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))[[1]]


icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  #filter(subdivisio != "20") %>% 
  as("Spatial") %>% 
  crop(extent(r))

m <- map_data("worldHires") 

m_sf <-  st_as_sf(x = m, 
                  coords = c("long", "lat"),
                  crs = crs(icesarea)) %>% 
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

m_sf <- as(m_sf, 'Spatial') %>% 
  crop(extent(Q1_layers_cut))

sob_c <- crop(mask(sob, icesarea), extent(r)) 

# Match plot specs to trawl coverage plot
png("../Figures/sob_fig.png", width = 17, height = 15, units = "cm", res = 300)
  # Bottom Salinity
plot(sob_c, col = oce.colorsSalinity(30), zlim = c(0,36),
     box=FALSE, axes =F, legend.width=1.5, axis.args=list(cex.axis=1.75), colNA = "grey70")
plot(icesarea, add = T, lwd = 1.25, interpolate = T)
#plot(m_sf, add = T, col="grey75", lwd = 0.75)
dev.off()

