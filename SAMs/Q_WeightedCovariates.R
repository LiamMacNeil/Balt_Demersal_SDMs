library(raster)
library(terra)
library(tidyverse)
library(oce)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read_csv("../Data/Cod_Herring_CatchWgt_Density.csv")

Quarter <- c("Q1", "Q4")
names(Quarter) <- c(1, 4)

# plot obs. frequencies for weighting
dat$Monthname <- lubridate::month(dat$Month, label = T)
dat %>%
  ggplot(aes(x= Monthname)) + 
  geom_bar()+  
  facet_wrap(~Quarter, scales = "free_x",
             labeller = labeller(Quarter=Quarter))+
  theme_bw(24)

# Get monthly weight vectors for each quarter

#Q1
dat %>%
  filter(Quarter ==1) %>% 
  group_by(Monthname) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# Feb = 0.446
# Mar = 0.554

Winter <- c(0.446, 0.554)

# Q4
dat %>%
  filter(Quarter ==4) %>% 
  group_by(Monthname) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# Oct = 0.181
# Nov = 0.806
# Dec 0.0136

Autumn <- c(0.181, 0.806, 0.0136)

################### 
################### Q1 Covariates
################### 

# Bottom Temp
bottomT <- stack(list.files("../../Predictors/Monthly/bottomT", 
                               pattern = "\\.tif$", full.names = TRUE))
bottomTQ1 <- bottomT[[c("February_bottomT", "March_bottomT")]]

bottomTQ1_mean <- mean(bottomTQ1)
bottomTQ1_w_mean <- weighted.mean(bottomTQ1, w = Winter)

#Chlorophyll
chl <- stack(list.files("../../Predictors/Monthly/chl", 
                            pattern = "\\.tif$", full.names = TRUE))
chlQ1 <- chl[[c("February_chl", "March_chl")]]

chlQ1_mean <- mean(chlQ1)
chlQ1_w_mean <- weighted.mean(chlQ1, w = Winter)

#MLD
mlotst <- stack(list.files("../../Predictors/Monthly/mlotst", 
                            pattern = "\\.tif$", full.names = TRUE))
mlotstQ1 <- mlotst[[c("February_mlotst", "March_mlotst")]]

mlotstQ1_mean <- mean(mlotstQ1)
mlotstQ1_w_mean <- weighted.mean(mlotstQ1, w = Winter)

#Surface O2
o2 <- stack(list.files("../../Predictors/Monthly/o2", 
                            pattern = "\\.tif$", full.names = TRUE))
o2Q1 <- o2[[c("February_o2", "March_o2")]]

o2Q1_mean <- mean(o2Q1)
o2Q1_w_mean <- weighted.mean(o2Q1, w = Winter)

# Bottom O2
o2b <- stack(list.files("../../Predictors/Monthly/o2b", 
                            pattern = "\\.tif$", full.names = TRUE))
o2bQ1 <- o2b[[c("February_o2b", "March_o2b")]]

o2bQ1_mean <- mean(o2bQ1)
o2bQ1_w_mean <- weighted.mean(o2bQ1, w = Winter)

# Surface Salinity
so <- stack(list.files("../../Predictors/Monthly/so", 
                            pattern = "\\.tif$", full.names = TRUE))
soQ1 <- so[[c("February_salin", "March_salin")]]

soQ1_mean <- mean(soQ1)
soQ1_w_mean <- weighted.mean(soQ1, w = Winter)

# Bottom Salinity
sob <- stack(list.files("../../Predictors/Monthly/sob", 
                            pattern = "\\.tif$", full.names = TRUE))
sobQ1 <- sob[[c("February_sob", "March_sob")]]

sobQ1_mean <- mean(sobQ1)
sobQ1_w_mean <- weighted.mean(sobQ1, w = Winter)

# Surface Temperature
thetao <- stack(list.files("../../Predictors/Monthly/thetao", 
                        pattern = "\\.tif$", full.names = TRUE))
thetaoQ1 <- thetao[[c("February_temp", "March_temp")]]

thetaoQ1_mean <- mean(thetaoQ1)
thetaoQ1_w_mean <- weighted.mean(thetaoQ1, w = Winter)

habitat <- raster("../../Predictors/CMEMS_Habitat_cut/GRIDCODE.tif")

# Reproject to combine
# Means
thetaoQ1_mean_reproj <- projectRaster(thetaoQ1_mean, habitat, method = "bilinear")
bottomTQ1_mean_reproj <- projectRaster(bottomTQ1_mean, habitat, method = "bilinear")
soQ1_mean_reproj <- projectRaster(soQ1_mean, habitat, method = "bilinear")
sobQ1_mean_reproj <- projectRaster(sobQ1_mean, habitat, method = "bilinear")
o2Q1_mean_reproj <- projectRaster(o2Q1_mean, habitat, method = "bilinear")
o2bQ1_mean_reproj <- projectRaster(o2bQ1_mean, habitat, method = "bilinear")
mlotstQ1_mean_reproj <- projectRaster(mlotstQ1_mean, habitat, method = "bilinear")
chlQ1_mean_reproj <- projectRaster(chlQ1_mean, habitat, method = "bilinear")

# Weighted Means
thetaoQ1_w_mean_reproj <- projectRaster(thetaoQ1_w_mean, habitat, method = "bilinear")
bottomTQ1_w_mean_reproj <- projectRaster(bottomTQ1_w_mean, habitat, method = "bilinear")
soQ1_w_mean_reproj <- projectRaster(soQ1_w_mean, habitat, method = "bilinear")
sobQ1_w_mean_reproj <- projectRaster(sobQ1_w_mean, habitat, method = "bilinear")
o2Q1_w_mean_reproj <- projectRaster(o2Q1_w_mean, habitat, method = "bilinear")
o2bQ1_w_mean_reproj <- projectRaster(o2bQ1_w_mean, habitat, method = "bilinear")
mlotstQ1_w_mean_reproj <- projectRaster(mlotstQ1_w_mean, habitat, method = "bilinear")
chlQ1_w_mean_reproj <- projectRaster(chlQ1_w_mean, habitat, method = "bilinear")

# Combine
 Q1_Mean_predictors <- stack(o2bQ1_mean_reproj, sobQ1_mean_reproj,
                           bottomTQ1_mean_reproj,chlQ1_mean_reproj,
                          habitat, mlotstQ1_mean_reproj,
                          o2Q1_mean_reproj, soQ1_mean_reproj, 
                          thetaoQ1_mean_reproj )
 names(Q1_Mean_predictors) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                        "Seabed Habitat", "MLD" ,"Surface.oxygen" , "Surface.salinity","Surface.temperature")
 Q1_w_Mean_predictors <- stack(o2bQ1_w_mean_reproj, sobQ1_w_mean_reproj,
                             bottomTQ1_w_mean_reproj,chlQ1_w_mean_reproj,
                             habitat, mlotstQ1_w_mean_reproj,
                             o2Q1_w_mean_reproj, soQ1_w_mean_reproj, 
                             thetaoQ1_w_mean_reproj )
 names(Q1_w_Mean_predictors) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                                "Seabed Habitat", "MLD" ,"Surface.oxygen" , "Surface.salinity","Surface.temperature")
 
 # Crop and resample to region excluding Gulf of Finland and Riga
 plot(crop(Q1_Mean_predictors, extent(8.987843, 22.5,
                                   53.91644, 59.25)))
 
 Q1_Mean_predictors <- crop(Q1_Mean_predictors, extent(8.987843, 22.5,
                                                       53.91644, 59.25))
 Q1_w_Mean_predictors <- crop(Q1_w_Mean_predictors, extent(8.987843, 22.5,
                                                       53.91644, 59.25))
 
 # Save Q1 mean stack
 for (p in 1:nlayers(Q1_Mean_predictors)) {
   writeRaster(Q1_Mean_predictors[[p]], filename = paste0("../../Predictors/Q1/",names(Q1_Mean_predictors)[p], ".tif"))
 }
 for (p in 1:nlayers(Q1_w_Mean_predictors)) {
   writeRaster(Q1_w_Mean_predictors[[p]], filename = paste0("../../Predictors/Q1/Weighted/",names(Q1_w_Mean_predictors)[p], ".tif"))
 }
 ################### 
 ################### Q4 Covariates
 ################### 
 
 # Bottom Temp
bottomTQ4 <- bottomT[[c("October_bottomT", "November_bottomT", "December_bottomT")]]
 
bottomTQ4_mean <- mean(bottomTQ4)
bottomTQ4_w_mean <- weighted.mean(bottomTQ4, w = Autumn)
 
#Chlorophyll
chlQ4 <- chl[[c("October_chl", "November_chl", "December_chl")]]
 
chlQ4_mean <- mean(chlQ4)
chlQ4_w_mean <- weighted.mean(chlQ4, w = Autumn)
 
 #MLD
mlotstQ4 <- mlotst[[c("October_mlotst", "November_mlotst", "December_mlotst")]]

mlotstQ4_mean <- mean(mlotstQ4)
mlotstQ4_w_mean <- weighted.mean(mlotstQ4, w = Autumn)

 #Surface O2
o2Q4 <- o2[[c("October_o2", "November_o2", "December_o2")]]

o2Q4_mean <- mean(o2Q4)
o2Q4_w_mean <- weighted.mean(o2Q4, w = Autumn)

 # Bottom O2
o2bQ4 <- o2b[[c("October_o2b", "November_o2b", "December_o2b")]]

o2bQ4_mean <- mean(o2bQ4)
o2bQ4_w_mean <- weighted.mean(o2bQ4, w = Autumn)

 # Surface Salinity
soQ4 <- so[[c("October_salin", "November_salin", "December_salin")]]

soQ4_mean <- mean(soQ4)
soQ4_w_mean <- weighted.mean(soQ4, w = Autumn)

 # Bottom Salinity
sobQ4 <- sob[[c("October_sob", "November_sob", "December_sob")]]

sobQ4_mean <- mean(sobQ4)
sobQ4_w_mean <- weighted.mean(sobQ4, w = Autumn)

 # Surface Temperature
thetaoQ4 <- thetao[[c("October_temp", "November_temp", "December_temp")]]

thetaoQ4_mean <- mean(thetaoQ4)
thetaoQ4_w_mean <- weighted.mean(thetaoQ4, w = Autumn)

 # Reproject to combine
 
# Means
 thetaoQ4_mean_reproj <- projectRaster(thetaoQ4_mean, habitat, method = "bilinear")
 bottomTQ4_mean_reproj <- projectRaster(bottomTQ4_mean, habitat, method = "bilinear")
 soQ4_mean_reproj <- projectRaster(soQ4_mean, habitat, method = "bilinear")
 sobQ4_mean_reproj <- projectRaster(sobQ4_mean, habitat, method = "bilinear")
 o2Q4_mean_reproj <- projectRaster(o2Q4_mean, habitat, method = "bilinear")
 o2bQ4_mean_reproj <- projectRaster(o2bQ4_mean, habitat, method = "bilinear")
 mlotstQ4_mean_reproj <- projectRaster(mlotstQ4_mean, habitat, method = "bilinear")
 chlQ4_mean_reproj <- projectRaster(chlQ4_mean, habitat, method = "bilinear")
 
 # Weighted Means
 thetaoQ4_w_mean_reproj <- projectRaster(thetaoQ4_w_mean, habitat, method = "bilinear")
 bottomTQ4_w_mean_reproj <- projectRaster(bottomTQ4_w_mean, habitat, method = "bilinear")
 soQ4_w_mean_reproj <- projectRaster(soQ4_w_mean, habitat, method = "bilinear")
 sobQ4_w_mean_reproj <- projectRaster(sobQ4_w_mean, habitat, method = "bilinear")
 o2Q4_w_mean_reproj <- projectRaster(o2Q4_w_mean, habitat, method = "bilinear")
 o2bQ4_w_mean_reproj <- projectRaster(o2bQ4_w_mean, habitat, method = "bilinear")
 mlotstQ4_w_mean_reproj <- projectRaster(mlotstQ4_w_mean, habitat, method = "bilinear")
 chlQ4_w_mean_reproj <- projectRaster(chlQ4_w_mean, habitat, method = "bilinear")
 
 # Combine
 Q4_Mean_predictors <- stack(o2bQ4_mean_reproj, sobQ4_mean_reproj,
                             bottomTQ4_mean_reproj,chlQ4_mean_reproj,
                             habitat, mlotstQ4_mean_reproj,
                             o2Q4_mean_reproj, soQ4_mean_reproj, 
                             thetaoQ4_mean_reproj )
 names(Q4_Mean_predictors) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                                "Seabed Habitat", "MLD" ,"Surface.oxygen" , "Surface.salinity","Surface.temperature")
 
 Q4_w_Mean_predictors <- stack(o2bQ4_w_mean_reproj, sobQ4_w_mean_reproj,
                             bottomTQ4_w_mean_reproj,chlQ4_w_mean_reproj,
                             habitat, mlotstQ4_w_mean_reproj,
                             o2Q4_w_mean_reproj, soQ4_w_mean_reproj, 
                             thetaoQ4_w_mean_reproj )
 names(Q4_w_Mean_predictors) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                                "Seabed Habitat", "MLD" ,"Surface.oxygen" , "Surface.salinity","Surface.temperature")
 
 # Crop and resample to region excluding Gulf of Finland and Riga
 plot(crop(Q4_Mean_predictors, extent(8.987843, 22.5,
                                      53.91644, 59.15)))
 
 Q4_Mean_predictors <- crop(Q4_Mean_predictors, extent(8.987843, 22.5,
                                                       53.91644, 59.25))
 Q4_w_Mean_predictors <- crop(Q4_w_Mean_predictors, extent(8.987843, 22.5,
                                                       53.91644, 59.25))
 # Save Q4 mean stack
 for (p in 1:nlayers(Q4_Mean_predictors)) {
   writeRaster(Q4_Mean_predictors[[p]], filename = paste0("../../Predictors/Q4/mean",names(Q4_Mean_predictors)[p], ".tif"))
 }
 for (p in 1:nlayers(Q4_w_Mean_predictors)) {
   writeRaster(Q4_w_Mean_predictors[[p]], filename = paste0("../../Predictors/Q4/Weighted/",names(Q4_w_Mean_predictors)[p], ".tif"))
 }
 
 # Extra (crude) plotting script to reproduce 
 plot(Q1_w_Mean_predictors[[1]], col = oce.colorsOxygen(20), zlim = c(-4,15))
 plot(Q4_w_Mean_predictors[[1]], col = oce.colorsOxygen(20), zlim = c(-4,12))
 
 plot(Q1_w_Mean_predictors[[7]], col = oce.colorsOxygen(20), zlim = c(5,15))
 plot(Q4_w_Mean_predictors[[7]], col = oce.colorsOxygen(20), zlim = c(5,12))
 
 plot(Q1_w_Mean_predictors[[9]], col = oce.colorsTemperature(20), zlim = c(0,12))
 plot(Q4_w_Mean_predictors[[9]], col = oce.colorsTemperature(20), zlim = c(0,12))
 
 plot(Q1_w_Mean_predictors[[3]], col = oce.colorsTemperature(20), zlim = c(0,9))
 plot(Q4_w_Mean_predictors[[3]], col = oce.colorsTemperature(20), zlim = c(3,14))
 
 plot(Q1_w_Mean_predictors[[8]], col = oce.colorsSalinity(20), zlim = c(5,34.5))
 plot(Q4_w_Mean_predictors[[8]], col = oce.colorsSalinity(20), zlim = c(5,34.5))
 
 plot(Q1_w_Mean_predictors[[2]], col = oce.colorsSalinity(20), zlim = c(5,36))
 plot(Q4_w_Mean_predictors[[2]], col = oce.colorsSalinity(20), zlim = c(5,36))
 
 plot(Q1_w_Mean_predictors[[4]], col = oce.colorsChlorophyll(20), zlim = c(0,4))
 plot(Q4_w_Mean_predictors[[4]], col = oce.colorsChlorophyll(20), zlim = c(0,4))
 
 plot(Q1_w_Mean_predictors[[6]], col = oce.colorsDensity(20))
 plot(Q4_w_Mean_predictors[[6]], col = oce.colorsDensity(20))
 