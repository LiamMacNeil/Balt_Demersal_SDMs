library(devtools)
library(tidyverse)
library(terra)
library(raster)
library(dismo)
library(oce)
library(gisland)
library(mapdata)
# install_version("spatstat", version = "1.64.1", repos = "http://cran.us.r-project.org")
library(spatstat) # need older version (1.64.1) for mopa
#library(mopa)
library(RColorBrewer)
library(viridis)
library(sjPlot)
library(gridExtra)
library(ggpubr)
library(blockCV)
library(sf)
library(mgcv)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################## ###################### ####################### #######################     
# Loading Data
####################### ####################### ####################### #######################

# Quarterly, Decadal Predictors
Q1_layers_cut <- stack(list.files("../Predictors/Monthly/2001-2010/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))

names(Q1_layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                          "MLD" ,"Seabed.Habitat","Surface.oxygen" , "Surface.salinity","Surface.temperature")

# Variable selection excluding surf. salinity + surf. oxygen displaying high multi-collinearity and VIF
varselect <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
               "MLD" ,"Seabed.Habitat","Surface.temperature")

# Pull boundary layers (https://github.com/einarhjorleifsson/gisland)
icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  filter(subdivisio != "20") %>% 
  as("Spatial") %>% 
  crop(extent(Q1_layers_cut))

m <- map_data("worldHires") 

m_sf <-  st_as_sf(x = m, 
                  coords = c("long", "lat"),
                  crs = crs(icesarea)) %>% 
  group_by(group) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 

m_sf <- as(m_sf, 'Spatial') %>% 
  crop(extent(Q1_layers_cut))

# Final crop of raster excluding Skagerrak
Q1_layers_balt <- crop(mask(Q1_layers_cut, icesarea), 
                       extent(icesarea))

# Raw quantities
dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv")

# Q1
dat <- dat %>% 
  #filter(ScientificName_WoRMS == "Pleuronectes platessa") %>% 
  #filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  #filter(ScientificName_WoRMS == "Gadus morhua") %>% 
  filter(ScientificName_WoRMS == "Limanda limanda") %>% 
  #filter(Year != 2008) %>%
  filter(HaulVal == "V" | HaulVal == "C"| HaulVal == "A") %>% 
  filter(Quarter == 1) %>% 
  filter(meanLong != "NA") %>% 
  filter(meanLat != "NA") %>%
  filter(Year < 2011) #%>%

# Size distinction used for Cod
# Adult
#filter(LngtClass > 349) 
# Juvenile
#filter(LngtClass < 350) 

#unique(dat$ScientificName_WoRMS)

# 2001-2010
#  Cod = 58,979
#  Flounder = 32,001
#  Plaice = 13,138
#  Dab = 10,512


# 2011-2020
#  Cod = 56,568
#  Flounder = 44,787
#  Plaice = 27,015
#  Dab = 13,142

######################## ###################### ####################### #######################     
# Creating predictor-response df 
####################### ####################### ####################### #######################

spatDat <- SpatialPointsDataFrame(dat[,c("meanLong","meanLat")], data.frame(dat[,"Density_g"]))
rast <- raster::rasterize(spatDat, 
                          Q1_layers_balt[[1]], 
                          field="Density_g")
pts <- rasterToPoints(rast,
                      fun = function(x){
                        x > 0
                      })

presence_cov <- raster::extract(Q1_layers_balt, pts[,1:2])

### Exlusion buffer method, reduces model fit for every species/time window
### Conceptually this is unsurprising given that sacrificing a greedy thinning procedure retains many more abundance observations which are highly variable
### Prevlaence is still accounted for under pseudo-absence generation
### Model predicitions are highly similar between methods, geographic distribution of "core abundance" areas is near identical but much lower when buffer applied
#bg <- backgroundGrid(Q1_layers_balt)
## Generate pseudo-absences considering an unique background extent
#RS_random <- pseudoAbsences(xy = data.frame(pts[,1:2]), 
#                            background = bg$xy, 
#                            exclusion.buffer = 0.15, 
#                            prevalence = 0.5, 
#                            kmeans = FALSE)

#absence <- (RS_random[["species1"]][["PA01"]][[1]])
#absence <- absence %>% 
#  dplyr::select(x,y)


#mask, no buffer
absence <- randomPoints(mask=Q1_layers_balt, n=nrow(pts),p=presence_cov)  

absence_cov <- raster::extract(Q1_layers_balt, absence)

presence_cov <- data.frame(presence_cov)
absence_cov <- data.frame(absence_cov)

presence_cov$Density_g <- pts[,3]
absence_cov$Density_g <- 0

presence_cov$Longitude <- pts[,1]
presence_cov$Latitude <- pts[,2]

absence_cov$Longitude <- absence[,1]
absence_cov$Latitude <- absence[,2]

full <- rbind(presence_cov, absence_cov)
full <- full[complete.cases(full),]
full <- full[!is.infinite(rowSums(full)),]

full$Density_kg <- (full$Density_g/1000 )
full$Density_log10 <- log10(full$Density_g + 1)

spc_col <- 14
var_cols <- 1:9  
names(full)[spc_col]  # check
names(full)[var_cols]  # check

full$Seabed.Habitat <- as.factor(full$Seabed.Habitat)
#######
# Data splitting
#######

# Spatial blocking CV
sf <- sf::st_as_sf(full, coords = c("Longitude", "Latitude"), crs = st_crs(icesarea))

# Spatial autocorrelarion of predictors
#png("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_blocking.png", width = 20, height = 12, res = 300, units = "cm")
sac1 <- blockCV::cv_spatial_autocor(r = Q1_layers_balt[[varselect]], 
                           num_sample = 5000)
#dev.off()

#png("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_variogram.png", width = 16, height = 10, res = 300, units = "cm")
plot(sac1$variograms[[1]])
#dev.off()

#ecv <- cv_cluster(x = test,
#                  column = "Density_log10",
#                  r = Q1_layers_balt[[varselect]],
#                  k = 5, 
#                  scale = TRUE)

#png("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_blocks.png", width = 16, height = 10, res = 300, units = "cm")
bloCV <- cv_spatial(x = sf,
                    column = "Density_log10",
                    size = sac1$range,
                    k = 5,
                    selection = "random",
                    iteration = 50)
#dev.off()

#png("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_MESS.png", width = 16, height = 10, res = 300, units = "cm")
cv_similarity(cv = bloCV, # the environmental clustering
              x = sf, 
              r = Q1_layers_balt[[varselect]], 
              progress = FALSE)
#dev.off()

#sb1 <- cv_spatial(x = test,
#                  column = "Density_log10", # the response column (binary or multi-class)
#                  k = 10, # number of folds
#                  size = sac1$range, # size of the blocks in metres
#                  selection = "random", # random blocks-to-fold
#                  iteration = 50) # also create folds for biomod2

#png("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_cv.png", width = 18, height = 12, res = 300, units = "cm")
cvp <- cv_plot(cv = bloCV, 
               x = sf,
               points_alpha = 0.5,
               label_size = 3)

#ggsave("../Figures/GAM/Spatial_autocorrelation/CodAdult_Q1_2020_cv.png", cvp, width = 24, height = 14, dpi = 300, units = "cm")

# Pulling full dataset with fold metadata
full_bCV <- cvp$data

#train <- sample(nrow(full), 0.75*nrow(full))
#dat_train <- full[train , ]
#dat_test <- full[-train, ]  
formula_GAM <- as.formula(paste(names(full_bCV)[12], "~", paste(varselect, collapse = "+")))
formula_GAM

gam1 <- mgcv::gam(Density_kg ~ Seabed.Habitat + s(Bottom.salinity)+s(Bottom.temperature)+s(Bottom.oxygen)+s(Chlorophyll)+s(MLD)+s(Surface.temperature), 
      data = full, family = gaussian(), method = 'ML', select = T, gamma = 1.4)
gam2 <- mgcv::gam(Density_kg ~ Seabed.Habitat + s(Bottom.salinity)+s(Bottom.temperature)+s(Bottom.oxygen)+s(Chlorophyll)+s(MLD)+s(Surface.temperature), 
                  data = full, family = tw(link = "log"), method = 'ML', select = T, gamma = 1.4)

AIC(gam1,gam2)

mae <- list()
cors <- list()
disp <- list()

folds <- unique(full_bCV$folds)

for (i in folds){
  
  dat_train_fold <- full_bCV %>% 
    filter(folds == i) %>% 
    filter(value == "Train")
  
  dat_test_fold <-  full_bCV %>%  
    filter(folds == i) %>% 
    filter(value == "Test")
  #dat_test_fold <- full %>% filter(fold_ID == i)
  
  #mod_GAM <- GAM(formula = formula_GAM, family = gaussian(link = "identity"), data = dat_train_fold)
  
  mod_GAM <- mgcv::gam(Density_kg ~ Seabed.Habitat + s(Bottom.salinity)+s(Bottom.temperature)+s(Bottom.oxygen)+s(Chlorophyll)+s(MLD)+s(Surface.temperature), 
                       data = full, family = tw(link = "log"), method = 'ML', select = T, gamma = 1.4)
  
  print(summary(mod_GAM))
  par(mfrow = c(2, 2))
  gam.check(mod_GAM)

  pred <- predict(mod_GAM, dat_test_fold)
  
  # MAE (accuracy)
  #print(mean(abs(dat_test_fold$Density_log10 - pred)))
  mae <- append(mae, list(mean(abs(dat_test_fold$Density_log10 - pred))))
  
  preds <- data.frame(dat_test_fold$folds,dat_test_fold$Density_log10, pred)
  
  # Discrimination (Pearson correlation - scaled and unscaled)
  cors <- append(cors, list(preds))
  
  # dispersion (precision) based on Waldock et al. 2022
  disp <- append(disp, list(  var(pred)/var(dat_test_fold$Density_log10)))
}

mae_df <- do.call(rbind, mae)
cors_df <- do.call(rbind, cors)
disp_df <- do.call(rbind, disp)

colnames(cors_df) <- c("folds","Observed_Density", "Predicted")

pred_cor <- cors_df %>% 
  #filter(Observed_Density > 0) %>% 
  #filter(Predicted > 0) %>% 
  ggplot(aes(x = (Observed_Density), y= Predicted))+
  geom_point(aes(color = folds), size = 2, alpha = 0.5)+
  geom_abline(slope=1, intercept=0)+
  theme_bw(16)+
  labs(x = "Log(Observed Density)", y = "Log(Predicted)")+
  stat_smooth(aes(color = folds), method=lm)+
  stat_cor(aes(color = folds, label = ..r.label..),
           method = "pearson", 
           size=4,  r.accuracy = 0.01, label.y.npc="top", label.x.npc = "left")
#guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)
ggsave("../Figures/GAM/Q1/CodAdult_Q1_2020_PredCor.png", pred_cor, width = 16, height = 12, units = "cm", dpi = 300)


discrimination <- list()

for (i in folds){
  
  cors_df_fold <- cors_df %>% filter(folds == i)
  
  covary <- cor(cors_df_fold$Observed_Density, cors_df_fold$Predicted, method = "pearson")
  
  discrimination <- append(discrimination, list(covary))
  
}

metrics <- data.frame(do.call(rbind, discrimination))
colnames(metrics) <- "discrimination"
metrics[,"mae"] <- data.frame(mae_df)
metrics[,"dispersion"] <- data.frame(disp_df)
metrics[,"folds"] <- data.frame(folds)

metrics$Quarter <- unique(dat$Quarter)
metrics$Year <- "2001-2010"
metrics$Taxa <- unique(dat$ScientificName_WoRMS)
metrics$model <- "GAM"

write.csv(metrics, "../Model_Metrics/GAM/CodAdult_Q1_2020_GAM.csv")

# Denisty plot (log transformed and normalized)
#pred_dens <- preds %>% 
#  filter(Observed_Density > 0) %>% 
#  filter(Predicted > 0) %>% 
#  ggplot(aes(x = ((Observed_Density-min(Observed_Density))/(max(Observed_Density)-min(Observed_Density))), 
#             y= ((Predicted-min(Predicted))/(max(Predicted)-min(Predicted))))) +
#  geom_density_2d_filled(alpha = 0.9, lwd = 0.25, color = 'black', bins=15)+
#  geom_point()+
#  geom_abline(slope = 1, intercept = 0, lwd = 2, linetype = 5)+
#  theme_bw(18)+
#  theme(legend.position = "none")+
#  labs(x = "Observed", y = "Predicted")+
#  scale_fill_viridis_d(option="viridis")
#ggsave("../Figures/GAM/Q1/CodAdult_Q1_2020_PredDens.png", pred_dens, width = 15, height = 12, units = "cm", dpi = 300)

# Deviance of residuals (Precision?)
1 - pchisq(deviance(mod_GAM), df.residual(mod_GAM))


############ Prediction ##################

GAM_final_mod <-mgcv::gam(Density_kg ~ Seabed.Habitat + s(Bottom.salinity)+s(Bottom.temperature)+s(Bottom.oxygen)+s(Chlorophyll)+s(MLD)+s(Surface.temperature), 
                          data = full, family = gaussian(), method = 'ML', select = T, gamma = 1.4)

par(mfrow = c(2, 2))
gam.check(GAM_final_mod)

#saveRDS(mod_GAM, file = paste0("../Models/GAM/CodAdult_Q1_2020_Valid.rds"))
saveRDS(GAM_final_mod, file = paste0("../Models/GAM/CodAdult_Q1_2020_Prediction.rds"))

par(mfrow=c(3, 3), mai = c(0.35, 0.35, 0.35, 0.35))

p <- plot_model(GAM_final_mod, type = "pred", title  = "")
nCol <- floor(sqrt((length(p))))

png("../Figures/GAM/Q1/CodAdult_Q1_2020_pdp.png", width = 15, height = 12, res = 300, units = "cm")
do.call("grid.arrange", c(p, ncol=nCol))
dev.off()

GAM_P <- predict(Q1_layers_balt[[varselect]], GAM_final_mod, type = "response")
GAM_P<- crop(mask(GAM_P, icesarea), extent(icesarea)) 

# Back-transformation (log10) to kg km^-2
GAM_P <- (10^(GAM_P)-1)/1000

GAM_P <- raster("../Predictions/CodAdult_Q1_2020_GAM.tif")

png("../Figures/GAM/Q1/CodAdult_Q1_2020_predmap.png", width = 15, height = 12, res = 300, units = "cm")
plot(abs(GAM_P),
     col = viridis(20), 
     axis.args=list(cex.axis=1.25),
     legend.width=1.5,
     
     zlim=c(0,0.5)
     
     #zlim=c(0,0.25)
     
     # Dab
     #zlim=c(0,0.25)
     
     # Plaice
     #zlim=c(0,0.075)
     
)
plot(icesarea, add = T, lwd = 1.25)
plot(m_sf, add = T, col="grey75", lwd = 0.75)
raster::contour(abs((GAM_P) > quantile(abs(GAM_P), probs = 0.95)), lwd = 0.25, col="firebrick2", drawlabels = F, add=T)
#points(dat_test %>% filter(Density_log10 >0) %>% dplyr::select(Longitude, Latitude) , pch=21, bg="purple")

dev.off()


writeRaster(GAM_P, "../Predictions/CodAdult_Q1_2020_GAM.tif", overwrite=TRUE)

