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
library(mopa)
library(fuzzySim)
library(RColorBrewer)
library(viridis)
library(rasterVis)
library(stars)
library(sf)
library(gbm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Quarterly Predictors
# 2001-2010
Q1_layers_cut <- stack(list.files("../Predictors/Monthly/2001-2010/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))


names(Q1_layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                          "MLD" ,"Seabed.Habitat","Surface.oxygen" , "Surface.salinity","Surface.temperature")

varselect <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
               "MLD" ,"Seabed.Habitat","Surface.oxygen" ,"Surface.temperature")

# Quantities
#dat <- read_csv("../Data/HH_HL_CA_Cod_Herring_Swept.csv")
dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv")

#Q1
dat <- dat %>% 
  filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  #filter(ScientificName_WoRMS == "Gadus morhua") %>% 
  #filter(Year != 2008) %>%
  filter(Quarter == 1) %>% 
  filter(meanLong != "NA") %>% 
  filter(meanLat != "NA") %>%
  filter(Year < 2011) 

#%>%
# Size distinction used for Cod
filter(LngtClass < 350) 

#  Cod
#  Q1 3935
#  Q4 3255

# Herring
#  Q1 3935
#  Q4 3412

spatDat <- SpatialPointsDataFrame(dat[,c("meanLong","meanLat")], data.frame(dat[,"Density_g"]))
rast <- raster::rasterize(spatDat, Q1_layers_cut[[1]], field="Density_g")
pts <- rasterToPoints(rast,
                      fun = function(x){
                        x > 0
                      })

presence_cov <- raster::extract(Q1_layers_cut, pts[,1:2])

absence <- randomPoints(mask=Q1_layers_cut, n=nrow(pts),p=presence_cov)  
absence_cov <- raster::extract(Q1_layers_cut, absence)

# 1178 equal points for each Q1 
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

full$Density_kg <- full$Density_g/1000 
full$Density_log10 <- log10(full$Density_g + 1)

spc_col <- 14
var_cols <- 1:9  
names(full)[spc_col]  # check
names(full)[var_cols]  # check

########## Grid search to optimize hyperparameters for retraining gbm.step to predict

train <- sample(nrow(full), 0.85*nrow(full))

dat_train <- full[train , ]
dat_test <- full[-train, ]  

# create hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(.001, .01, .1),
  interaction.depth = c(1, 5, 8),
  n.minobsinnode = c(5, 10, 15),
  bag.fraction = c(.25, .5, .75), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# total number of combinations
nrow(hyper_grid)

# ranbdomize data order first:
random_index <- sample(1:nrow(dat_train), nrow(dat_train))
dat_train <- dat_train[random_index, ]
formula_GBM <- as.formula(paste(names(dat_train)[spc_col], "~", paste(varselect, collapse = "+")))

# Run grid search against parameters 
# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(123)
  
  # train model
  gbm.tune <- gbm(
    formula = formula_GBM,
    distribution = "tweedie",
    data = dat_train,
    n.trees = 5000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .8,
    cv.folds = 10,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

# Find specs for optimized model params
hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(15) 

pred_tune <- predict(gbm.tune, dat_test, 321)

preds_tune <- data_frame(dat_test[,spc_col], pred_tune)
colnames(preds_tune) <- c("Observed_Density", "Predicted")

preds_tune %>% filter(Observed_Density > 0) %>% 
  filter(Predicted > 0) %>% 
  ggplot(aes(x = ((Observed_Density-min(Observed_Density))/(max(Observed_Density)-min(Observed_Density))), 
             y= ((Predicted-min(Predicted))/(max(Predicted)-min(Predicted))))) +
  geom_density_2d_filled(alpha = 0.9, lwd = 0.25, color = 'black', bins=15)+
  geom_point(color="white")+
  geom_abline(slope = 1, intercept = 0, lwd = 2, linetype = 5)+
  theme_bw(18)+
  theme(legend.position = "none")+
  labs(x = "Observed", y = "Predicted")+
  scale_fill_viridis_d(option="plasma")
