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

icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  as("Spatial")

# Quarterly Predictors
# 2001-2010
Q1_layers_cut <- stack(list.files("../Predictors/Monthly/2001-2010/Q1/Weighted/", 
                                  pattern = "\\.tif$", full.names = TRUE))


names(Q1_layers_cut) <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
                          "MLD" ,"Seabed.Habitat","Surface.oxygen" , "Surface.salinity","Surface.temperature")

varselect <- c("Bottom.oxygen" ,"Bottom.salinity", "Bottom.temperature","Chlorophyll" ,
               "MLD" ,"Seabed.Habitat","Surface.oxygen" ,"Surface.temperature")

# Quantities
#dat <- read_csv("../Data/Cod_Herring_CatchWgt_Density.csv")
#dat <- read_csv("../Data/Cod_Herring_CatchWgt_Density_doorswept.csv")
dat <- read_csv("../Data/HH_HL_CA_Cod_Herring_Swept.csv")
dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv")

#Q1
dat <- dat %>% 
  filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  #filter(ScientificName_WoRMS == "Gadus morhua") %>% 
  #filter(Year != 2008) %>%
  filter(Quarter == 1) %>% 
  filter(meanLong != "NA") %>% 
  filter(meanLat != "NA") %>%
  filter(Year < 2011) #%>%
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

########################  
####################### GENERALIZED LINEAR MODEL (GLM)
####################### 

formula_GLM <- as.formula(paste(names(full)[spc_col], "~", paste(varselect, collapse = "+")))
formula_GLM

mod_GLM <- glm(formula = formula_GLM, family = gaussian, data = full)
summary(mod_GLM)
full$GLM_P <- predict(mod_GLM, newdata = full, type = "response")

AIC(mod_GLM)
BIC(mod_GLM)

#Deviance of residuals 
1 - pchisq(deviance(mod_GLM), df.residual(mod_GLM))

####################### 
####################### BOOSTED REGRESSION TREES / GENERALIZED BOOSTED MODELS (GBM)
####################### 

############ Validation ##################
train <- sample(nrow(full), 0.8*nrow(full))

dat_train <- full[train , ]
dat_test <- full[-train, ]  

model_bf <- gbm.step(dat_train, gbm.x = varselect, gbm.y = spc_col,
                     family = "gaussian", 
                     tree.complexity = 10, 
                     learning.rate = 0.005, 
                     bag.fraction = 0.7, 
                     plot.folds = T)

# Plot model output
gbm.plot(model_bf)


simply <- gbm.simplify(model_bf,  plot = T)

model_bf_simp <- gbm.step(dat_train, 
                          gbm.x = simply$pred.list[[1]], 
                          gbm.y = spc_col,
                          family = "gaussian", 
                          tree.complexity = 5, 
                          learning.rate = 0.005, 
                          #n.folds = 10,
                          bag.fraction = 0.5#, plot.folds = T
                          )

gbm.plot(model_bf_simp)

# Plot fitted vs. observed values (wtm = weighted mean of fitted values in relation to each pred)
gbm.plot.fits(model_bf_simp)

# Check if pairwise interactions exist
find.int <- gbm.interactions(model_bf_simp)
find.int$interactions #mainly interactions with vars calculated from each other

# Plot intercations
gbm.perspec(model_bf_simp, 1, 7, z.range=c(0,2))

model_bf_simp$cv.statistics

pred <- predict(model_bf_simp, dat_test, model_bf_simp$gbm.call$best.trees)

cor(pred, dat_test[spc_col])

mean(abs(dat_test[,spc_col]-pred))

plot(dat_test[,spc_col], pred)

preds <- data_frame(dat_test[,spc_col], pred)
colnames(preds) <- c("Observed_Density", "Predicted")

preds %>% filter(Observed_Density > 0) %>% 
  filter(Predicted > 0) %>% 
  ggplot(aes(x = Observed_Density, y= Predicted))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_abline(slope=1, intercept=0)+
  theme_bw(18)+
  labs(x = "Log(Observed Density)", y = "Predicted")

# Denisty plot (log transformed and normalized)
pred_dens <- preds %>% filter(Observed_Density > 0) %>% 
  filter(Predicted > 0) %>% 
  ggplot(aes(x = ((Observed_Density-min(Observed_Density))/(max(Observed_Density)-min(Observed_Density))), 
             y= ((Predicted-min(Predicted))/(max(Predicted)-min(Predicted))))) +
  geom_density_2d_filled(alpha = 0.9, lwd = 0.25, color = 'black', bins=15)+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, lwd = 2, linetype = 5)+
  theme_bw(18)+
  theme(legend.position = "none")+
  labs(x = "Observed", y = "Predicted")+
  scale_fill_viridis_d(option="plasma")
ggsave("../Figures/Her_Q1Pred_dens.png", pred_dens, width = 15, height = 12, units = "cm", dpi = 300)


# Model's deviance
# source: https://stats.stackexchange.com/questions/108995/interpreting-residual-and-null-deviance-in-glm-r
D <- model_bf_simp$self.statistics$mean.resid #residual deviance
n <- length(model_bf_simp$var.levels) #number of regressors
p <- nrow(full) #number of observations
D/(n - p) #if >> 1, rule of thumb is that model is inadequate

saveRDS(model_bf_simp, file = paste0("../Models/Q1Cod_Full.rds"))

###########
### Consult sensitivity script (GBM_Sensitivity)
###########

############ Prediction ##################

model_bf_simp_final <- gbm.step(full, gbm.x=simply$pred.list[[1]], gbm.y = spc_col,
                          family = "gaussian",
                          tree.complexity = 5, 
                          learning.rate = 0.005, 
                          #n.folds = 10,
                          bag.fraction = 0.8, 
                          plot.folds = T)

GBM_P <- predict(Q1_layers_cut[[simply$pred.list[[1]]]], model_bf_simp_final, type = "response", 
                 n.trees=model_bf_simp_final$gbm.call$best.trees)
 
GBM_P<- crop(mask(GBM_P, icesarea), extent(icesarea)) 

png("../Figures/Q1Her_predmap.png", width = 15, height = 12, res = 300, units = "cm")
plot(abs(GBM_P),
     col = plasma(20), 
     axis.args=list(cex.axis=1.25),
     legend.width=1.5,
     zlim=c(0,3))
contour(abs((GBM_P) > quantile(abs(GBM_P), probs = 0.9)), lwd = 0.00001, drawlabels = F, add=T)
dev.off()

saveRDS(model_bf_simp_final, file = paste0("../Models/Q1Her_Prediction.rds"))

#################
################# Habitat Overlap
#################

# Cod
cod_spawn <- (raster("../../../Cod-Her-SDMs/Cod_Spawning_HELCOM/EFH_Cod.tif"))
cod_spawn <- projectRaster(cod_spawn,GBM_P,method = 'bilinear')

cod_spawn <- crop(mask(cod_spawn, GBM_P), extent(GBM_P))

plot(cod_spawn, legend = FALSE, col = rev(terrain.colors(3)), box=F, axes=F)
legend("topright", legend = c("NA","Probable", "Essential"), fill = rev(terrain.colors(3)), cex=2)

# Herring
herring_spawn <- (raster("../../../Cod-Her-SDMs/Herring_Spawning_HELCOM/EFH_Herring.tif"))
herring_spawn <- projectRaster(herring_spawn,GBM_P,method = 'bilinear')

herring_spawn <- crop(mask(herring_spawn, GBM_P), extent(GBM_P))

plot(herring_spawn, legend = FALSE, col = rev(terrain.colors(3)), box=F, axes=F)
legend("topright", legend = c("NA","Probable", "Essential"), fill = rev(terrain.colors(3)), cex=2)

plot(abs(GBM_P), legend=F, , xlab=  expression('Longitude '*degree*"E"),
     ylab=expression("Latitude "*degree*"N"), 
     cex.axis=1.3, cex.main=1.8, cex.lab=1.6)

# Overlap two rasters and reclassify NAs to not appear
masked <- mask(abs((GBM_P) > quantile(abs(GBM_P), probs = 0.9)), 
               reclassify(cod_spawn>0, cbind(-Inf, 0, NA), right=T))

masked <- reclassify(masked, cbind(-Inf, 0, NA), right=T)

plot(masked, col="dodgerblue1", add=T, legend=F)


plot(abs((GBM_P) > quantile(abs(GBM_P), probs = c(0.1, 0.9))[2]), 
     col = hcl.colors(20), legend = NULL)

# map data
xlim <- c(9.25, 26.5)
ylim <- c(53.5, 58.5)

m <- map_data("worldHires", xlim = xlim, ylim = ylim)

coords <- xyFromCell(GBM_P[[1]], seq_len(ncell(GBM_P[[1]])))
df <- stack(as.data.frame(getValues(GBM_P[[1]])))
names(df) <- c('value', 'variable')
df$variable <- "Predicted"
df <- cbind(coords, df)

df %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = value))+
  geom_polygon(data = m, aes(long, lat, group = group), fill = "white") +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  #geom_tile(aes(x, y, fill = value), data = df, inherit.aes = T)+
  geom_path(data = icesarea, aes(long, lat, group = group), size=1) +
  theme_bw(18) +
  xlab(expression('Longitude '*degree*"E"))+
  ylab(expression("Latitude "*degree*"N"))+
  scale_fill_gradientn(colours = plasma(10), name=expression(log(kg~km^2)))


############## 
############## Extra stuff (under development)
############## 

GBM_P_sf <- st_as_sf(st_as_stars(GBM_P))
m_sf <- st_as_sf(m, coords = c("long", "lat"), crs = st_crs(GBM_P_sf))

GBM_P_sf %>% 
  ggplot() +
  geom_sf(aes(fill=layer), lwd=0)+
  #geom_sf(data = m_sf, fill = "white") +
  coord_sf()+ 
  theme_bw(20)+
  xlab("")+
  ylab("")


GLM_P <- predict(Q1_layers_cut, mod_GLM, type = "response")
plot(abs(GLM_P), col = hcl.colors(20), main = "GLM Raw Density (2001-2020)")

#  Predictions on data table:
full$GBM_P <- predict(mod_GBM, newdata = full, type = "response", n.trees = 100)

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
formula_RF <- as.formula(paste("Density ~", paste(varselect, collapse = "+")))
mod_RF <- randomForest(formula = formula_RF, data = full, na.action = na.exclude)  # I've used the defaults, but read the help file to explore important parameters like 'ntree', 'mtry', etc. if you use RF for real work; see also http://uc-r.github.io/2018/05/09/random-forests/
mod_RF
# get predictions on Raster:
RF_P <- predict(Q1_layers_cut, mod_RF, type = "response")
plot(RF_P, col = hcl.colors(20), main = "RF")


?gbm
formula_GBM <- formula_GLM
mod_GBM <- gbm(formula_GBM, distribution = "gaussian" , data = full,
               n.trees = 5000, cv.folds = 10, interaction.depth = 5,
               shrinkage = 0.05) 

summary(mod_GBM)

# find index for n trees with minimum CV error
min_MSE <- which.min(mod_GBM$cv.error)
# get MSE and compute RMSE
sqrt(mod_GBM$cv.error[min_MSE])

# plot loss function as a result of n trees added to the ensemble
gbm.perf(mod_GBM, method = "cv", oobag.curve = T)

auc(mod_GBM)

