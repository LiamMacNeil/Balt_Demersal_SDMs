library(tidyverse)
library(terra)
library(raster)
library(dismo)
library(oce)
library(embarcadero)
library(gisland)
library(mapdata)
library(mopa)
library(fuzzySim)
library(RColorBrewer)
library(rasterVis)
library(gam)
library(dismo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

layers_cut <- stack(list.files("C:/Users/lmacneil/Documents/copernicus-processed-data/covariates/CMEMS_Habitat_cut", 
                               pattern = "\\.tif$", full.names = TRUE))

names(layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                       "Seabed Habitat", "MLD" ,"Surface.oxygen" , "Surface.salinity","Surface.temperature")


dat <- read_csv("C:/Users/lmacneil/Documents/DATRAS/StandardWgt_Valid_and_Pelagic.csv")

# Extreme haul duration values?
dat <- dat %>% 
  filter(HaulDur >! 90) %>% 
  filter(HaulDur >= 15)

# Calculating distance manually? 
dat <- dat %>% 
  filter(Distance !="NA") %>% 
  filter(Taxa == "Cod") %>% 
  filter(Year > 2001)

dat <- dat %>% 
  dplyr::mutate(perarea = ifelse(DataType=="R" | DataType == "S",
                                 CatchWgt / Distance,
                                 CatchWgt))


#dat <- dat[,c("ShootLong", "ShootLat",
#              "CatchWgt",
#              "Year")]

spatDat <- SpatialPointsDataFrame(dat[,c("ShootLong","ShootLat")], data.frame(dat[,"perarea"]))
min <- raster::rasterize(spatDat, layers_cut[[1]], field="perarea", fun="min")
pts_min <- rasterToPoints(min,
                          fun = function(x){
                            x > 0
                          })

presence_cov <- raster::extract(layers_cut, pts_min[,1:2])

absence <- randomPoints(mask=layers_cut, n=nrow(pts_min),p=presence_cov)  
absence_cov <- raster::extract(layers_cut, absence)

# 1987 equal points for each presence and absence 
presence_cov <- data.frame(presence_cov)
absence_cov <- data.frame(absence_cov)

presence_cov$Abundance <- pts_min[,3]
absence_cov$Abundance <- 0

presence_cov$Longitude <- pts_min[,1]
presence_cov$Latitude <- pts_min[,2]

absence_cov$Longitude <- absence[,1]
absence_cov$Latitude <- absence[,2]

full <- rbind(presence_cov, absence_cov)
full <- full[complete.cases(full),]
full <- full[!is.infinite(rowSums(full)),]

full$Abundance <- log10(full$Abundance +1)
  
spc_col <- 10
var_cols <- 1:9  
names(full)[spc_col]  # check
names(full)[var_cols]  # check


varselect <- c("Bottom.oxygen" ,"Bottom.salinity" ,"Bottom.temperature","Chlorophyll" ,"Seabed.Habitat",
               "MLD" ,"Surface.oxygen","Surface.temperature")

# GENERALIZED LINEAR MODEL (GLM)

formula_GLM <- as.formula(paste(names(full)[spc_col], "~", paste(varselect, collapse = "+")))
formula_GLM
mod_GLM <- glm(formula = formula_GLM, family = gaussian, data = full)
summary(mod_GLM)
full$GLM_P <- predict(mod_GLM, newdata = full, type = "response")


# BOOSTED REGRESSION TREES / GENERALIZED BOOSTED MODELS (GBM)
train <- sample(nrow(full), 0.7*nrow(full))

dat_train <- full[train , ]
dat_test <- full[-train, ]  

learningrates <- 0.003
bagfrac <- 0.6

model_bf <- gbm.step(full, gbm.x = varselect, gbm.y = spc_col,
                     family = "gaussian", tree.complexity = 3, learning.rate = learningrates,
                     bag.fraction = bagfrac)

pred <- predict(model_bf, dat_test, model_bf$gbm.call$best.trees)

cor(pred, dat_test[spc_col])

mean(abs(dat_test[,spc_col]-pred))

plot(dat_test[,spc_col], pred)

# Model's deviance
# source: https://stats.stackexchange.com/questions/108995/interpreting-residual-and-null-deviance-in-glm-r
D <- model_bf$self.statistics$mean.resid #residual deviance
n <- length(model_bf$var.levels) #number of regressors
p <- nrow(full) #number of observations
D/(n - p) #if >> 1, rule of thumb is that model is inadequate

# Plot model output
gbm.plot(model_bf, n.plots = n, plot.layout = c(3, 3), write.title = F)

# Plot fitted vs. observed values (wtm = weighted mean of fitted values in relation to each pred)
gbm.plot.fits(model_bf)


# Check if pairwise interactions exist
find.int <- gbm.interactions(model_bf)
find.int$interactions #mainly interactions with vars calculated from each other


?gbm
formula_GBM <- formula_GLM
mod_GBM <- gbm(formula_GBM, distribution = "gaussian" , data = full)  # I've used the defaults, but read the help to explore important parameters like 'shrinkage' (learning rate) and 'interaction.depth' (tree complexity) if you use GBM/BRT for real work; see also https://doi.org/10.1111/j.1365-2656.2008.01390.x
summary(mod_GBM)
# get predictions on the data table:
full$GBM_P <- predict(mod_GBM, newdata = full, type = "response", n.trees = 100)


# APPLY PREDICTIONS TO THE RASTER VARIABLES ####

GLM_P <- predict(layers_cut, mod_GLM, type = "response")
plot(abs(GLM_P), col = hcl.colors(20), main = "GLM Raw Abundance (2001-2021)")


GBM_P <- predict(layers_cut, mod_GBM, type = "response", n.trees = 100)
plot(abs(GBM_P), col = hcl.colors(20), main = "GBM Raw Abundance (2001-2021)")


pred_rasters <- stack(GLM_P, GAM_P, RF_P, GBM_P)

plot(pred_rasters,  col=hcl.colors(10),  cex = 2)

plot(mean(pred_rasters), col=hcl.colors(10))

GAM_P <- predict(layers_cut, mod_GAM, type = "response")
plot(GAM_P, col = hcl.colors(20), main = "GAM")

RF_P <- predict(layers_cut, mod_RF)
plot(RF_P, col = hcl.colors(10), main = "RF")
# GENERALIZED ADDITIVE MODEL (GAM)

?gam
formula_GAM <- as.formula(paste(names(full)[spc_col], "~", paste0("s(", varselect, ")", collapse = "+")))  # GAM with smoothing splines ('s')
formula_GAM
mod_GAM <- gam(formula = formula_GAM, family = gaussian, data = full)
summary(mod_GAM)
# get predictions on the data table:
full$GAM_P <- predict(mod_GAM, newdata = full, type = "response")

# RANDOM FOREST (RF)
?randomForest
formula_RF <- as.formula(paste("Abundance ~", paste(varselect, collapse = "+")))
mod_RF <- randomForest(formula = formula_RF, data = full, na.action = na.exclude)  # I've used the defaults, but read the help file to explore important parameters like 'ntree', 'mtry', etc. if you use RF for real work; see also http://uc-r.github.io/2018/05/09/random-forests/
mod_RF
# get predictions on the data table:
full$RF_P <- predict(mod_RF, newdata = full)

