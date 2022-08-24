library(gurobi)
library(raster)
library(rgdal)
library(prioritizr)
library(prioritizrdata)
library(oce)
library(scales)
library(withr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  as("Spatial")

Cod <- raster("../BO_Cod/Static/Outputs/Predicted_Static.tif")
Herring <- raster("../BO_Herring/Static/Outputs/Predicted_Static.tif")
features <- brick(Cod, Herring)
names(features) <- c("Cod", "Herring")

features_trim <- features>0.5

Velocity <- raster("VoCC/voccMag.tif")
Velocity <- crop(mask(Velocity, icesarea), extent(icesarea))
Velocitytresampled <- projectRaster(Velocity,features,method = 'bilinear')

MPAs <- raster("../Balt_MPAs/MPA_rasters/MPA_ID.tif")
MPArast_resamp <- crop(mask(MPAs, icesarea), extent(Velocitytresampled))
MPArast_resamp <- projectRaster(MPArast_resamp,features,method = 'bilinear')
MPArast_resamp <- MPArast_resamp >0

boundary_data <- boundary_matrix(Velocitytresampled)
#boundary_data@x <- rescale(boundary_data@x, to = c(0.01, 100))

# Baseline prioritization
p <- problem(Velocitytresampled, features_trim) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.10) %>% 
  add_boundary_penalties(penalty = 80, data = boundary_data) %>% 
  #add_locked_in_constraints(MPArast_resamp) %>% 
  add_binary_decisions() %>% 
  add_gurobi_solver() 

p1 <- problem(Velocitytresampled, features_trim) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.25) %>% 
  add_boundary_penalties(penalty = 80, data = boundary_data) %>% 
  #add_locked_in_constraints(MPArast_resamp) %>% 
  add_binary_decisions() %>% 
  add_gurobi_solver() 

p2 <- problem(Velocitytresampled, features_trim) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.5) %>% 
  add_boundary_penalties(penalty = 80, data = boundary_data) %>% 
  #add_locked_in_constraints(MPArast_resamp) %>% 
  add_binary_decisions() %>% 
  add_gurobi_solver() 

p3 <- problem(Velocitytresampled, features_trim) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.75) %>% 
  add_boundary_penalties(penalty = 80, data = boundary_data) %>% 
  #add_locked_in_constraints(MPArast_resamp) %>% 
  add_binary_decisions() %>% 
  add_gurobi_solver() 

s <- stack(lapply(list(p, p1,p2, p3), solve))

plot(s, main = c("10%","25%","50%", "75%"), 
     axes = T, legend=F,col = c("grey60", "forestgreen"), cex.axis=1.5, cex.main=1.5,
)

######### Useless as MPA locking drives 1:1 between prioritized and irreplaceable (seemingly only solution)

px <- stack(eval_replacement_importance(p, s$voccMag.1),
            eval_replacement_importance(p1, s$voccMag.2),
            eval_replacement_importance(p2, s$voccMag.3),
            eval_replacement_importance(p3, s$voccMag.4))

px <- crop(mask(px, icesarea), extent(icesarea))

# set infinite values as 1.09 so we can plot them
px$rc.1[px$rc.1 > 100] <- 1.09
px$rc.2[px$rc.2 > 100] <- 1.09
px$rc.3[px$rc.3 > 100] <- 1.09
px$rc.4[px$rc.4 > 100] <- 1.09

# plot the importance scores

plot(px)
# planning units that are truly irreplaceable are shown in red
par(mfrow = c(2,2))
plot(px, main=c("10% Target", "25% Target",
                "50% Target", "75% Target"),
     at = c(seq(0, 0.9, 0.1), 1.01, 1.1), col=oce.colorsPalette(10),
     cex.main=1.8)

#col = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
    #         "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725",
    #         "#FF0000"))
