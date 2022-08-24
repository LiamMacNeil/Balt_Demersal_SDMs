#install.packages("C:/gurobi952/win64/R/gurobi_9.5-2.zip", repos = NULL)
#install.packages("prioritizr", repos = "https://cran.rstudio.com/")
#install.packages("slam", repos = "https://cloud.r-project.org")

library(gurobi)
library(raster)
library(terra)
library(rgdal)
library(prioritizr)
library(prioritizrdata)
library(oce)
library(scales)
library(withr)
library(tidyverse)
library(ENMTools)

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

#grid <- rasterToPolygons(Velocitytresampled)
#climate_connect <- raster("VoCC/ClassL_Velocity.tif")
#climate_connect <- raster("VoCC/TrajClas_Velocity.tif")
#climate_connect <- crop(mask(climate_connect, icesarea), extent(icesarea))
#climate_connectresampled <- projectRaster(climate_connect,features,method = 'bilinear')

#Balt <- read_sf("../Balt_shape" , layer = "Baltic")
#windfarm <- read_sf("C:/Users/lmacneil/Documents/outputs/Balt_WindFarms" , layer = "Wind_turbines")
MPAs <- raster("../Balt_MPAs/MPA_rasters/MPA_ID.tif")
MPArast_resamp <- crop(mask(MPAs, icesarea), extent(Velocitytresampled))
MPArast_resamp <- projectRaster(MPArast_resamp,features,method = 'bilinear')
MPArast_resamp <- MPArast_resamp >0


#############
############# UNLOCKED CONSTRAINTS (NO MPAS CONSIDERED)
#############

# Baseline prioritization
p <- problem(Velocitytresampled, features_trim) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.25) %>% 
  add_binary_decisions() %>% 
  add_gurobi_solver() 
  
s <- solve(p)

plot(
  s, main = "Baseline prioritization", axes = FALSE, legend=F,
  breaks = c(0, 0.5, 1), col = c("grey60", "forestgreen")
)

#

boundary_data <- boundary_matrix(Velocitytresampled)
#boundary_data@x <- rescale(boundary_data@x, to = c(0.01, 100))

p1 <- p %>%
  add_connectivity_penalties(penalty = 20, data = connectivity_matrix(Velocitytresampled, features[[1]]))

p2 <- p %>%
  # Creates snaking, discontiuous areas
  #add_neighbor_constraints(k=2) %>% 
  # Extraordinarily complex, not implementable
  #add_feature_contiguity_constraints() %>% 
  #add_linear_constraints(threshold = threshold, sense = "<=", data = climate_connectresampled) %>%
  add_boundary_penalties(penalty = 80, data = boundary_data) 

p3 <- p2 %>% 
  add_locked_in_constraints(MPArast_resamp) 

#p3 <- p2 %>%  add_feature_contiguity_constraints()

s <- stack(lapply(list(p, p1,p2, p3), solve))

s <- crop(mask(s, icesarea), extent(icesarea))

plot(s, main = c("Baseline","Connectivity Penalty","Boundary Penalty", "Locked-in MPAs"), 
     axes = T, legend=F,col = c("grey60", "forestgreen"),
     cex.axis=1.5, cex.main=1.5
)

#Evaluate moundary length modifier (already specified in penalty)
#r1 <- eval_boundary_summary(p1, s1)

#############
############# Intersection 
#############

#from Warren et al. 2008 
raster.overlap(features_trim[[1]], features_trim[[2]])


plot(s$voccMag.4, legend=F, col=c("grey60","forestgreen"), xlab=  expression('Longitude '*degree*"E"),
       ylab=expression("Latitude "*degree*"N"), main= "10% Target", 
     cex.axis=1.3, cex.main=1.8, cex.lab=1.6)
plot(MPArast_resamp, legend=F, col="dodgerblue4", add=T)

# Overlap two rasters and reclassify NAs to not appear
masked <- mask(s$voccMag.3, MPArast_resamp)
masked <- reclassify(masked, cbind(-Inf, 0, NA), right=T)

plot(masked, col="dodgerblue1", add=T, legend=F)

#############
############# Sensitivity 
#############

# Layers will be converted to spatialpolygons for interoperability
Velocitytresampled_sf <-st_as_sf(as(Velocitytresampled, "SpatialPolygonsDataFrame"))
#features_sf <-st_as_sf(as(features, "SpatialPolygonsDataFrame"))
MPArast_resamp_sf <-st_as_sf(as(MPArast_resamp, "SpatialPolygonsDataFrame"))

Baltic_PU <- Velocitytresampled_sf
Baltic_PU$PU <- 0

Baltic_PU <- st_join(Baltic_PU, MPArast_resamp_sf)
names(Baltic_PU)<- c("voccMag" , "PU"     ,  "MPAs"   , "geometry")

Baltic_PU[is.na(Baltic_PU)] <- 0

###
# Blended Approach
lower <- 0.75
upper <- 2.5
penalty <- round(10^seq(lower, upper, length.out=9),1)

boundary_data <- boundary_matrix(Baltic_PU)
boundary_data@x <- rescale(boundary_data@x, to = c(0.01, 100))

Locked_in <- Baltic_PU[Baltic_PU$MPAs==1,]

pblend <- problem(Baltic_PU, features_trim, cost_column="voccMag") %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.5) %>% 
  #add_locked_in_constraints(Locked_in[,"MPAs"]) %>% 
  add_binary_decisions() 

blended <- lapply(penalty, function(x){
  sblend <-  pblend %>% 
    add_boundary_penalties(penalty = x, data = boundary_data) %>% 
    add_gurobi_solver()  %>% 
    solve()
  s <- as.data.frame(sblend$solution_1)
  names(sblend) <- with_options(list(scipen=30), paste0("penalty_",x))
  s
})

blended_results <- cbind(
  Baltic_PU, do.call(bind_cols, blended)
)

colnames(blended_results)[4:12] <-  c(paste0("Penalty ",penalty))

svg("Prioritization_50_BoundaryLength_Sensitivity.svg", height = 16, width = 22)
# plot maps of prioritizations
plot(
  x =
    blended_results %>%
    dplyr::select(starts_with("penalty ")) %>%
    mutate_if(is.numeric, function(x) {
      case_when(
        blended_results$MPAs > 0.5 ~ "locked in",
        x > 0.5 ~ "priority",
        TRUE ~ "other"
      )
    }),
  pal = c("lightslateblue", "white", "forestgreen"), cex.main = 2.15, lwd=0.1
)
dev.off()

blended_metrics <- lapply(
  grep("Penalty ", names(blended_results)), function(x) {
    x <- blended_results[ , x]
    data.frame(
      total_cost = eval_cost_summary(pblend , x)$cost,
      total_boundary_length = eval_boundary_summary(pblend , x)$boundary
    )
  }
)

blended_metrics <- do.call(bind_rows, blended_metrics)
blended_metrics$penalty <- penalty
blended_metrics <- as_tibble(blended_metrics)

result_data <-
  blended_metrics %>%
  ## rename threshold column to value column
  rename(value = "penalty") %>%
  ## add column with column names that contain candidate prioritizations
  mutate(name = grep(
    "penalty", names(blended_metrics), value = TRUE, fixed = TRUE
  )) %>%
  ## add column with labels for plotting
  mutate(label = paste(value)) %>%
  ## add column to keep track prioritizations selected by different methods
  mutate(method = "none")

result_plot <-
  ggplot(
    data = result_data,
    aes(x = total_boundary_length, y = total_cost)
  ) +
  geom_line(size=1) +
  geom_point(aes(colour=factor(label)) ,size = 4) +
  scale_color_viridis_d()+
  xlab("Total Boundary Length") +
  ylab("Total Cost (VoCC)") +
  theme_bw(20)+
  guides(color=guide_legend(title="Boundary Penalty"))

ggsave("Penalty_Cost_BoundaryLength.svg" , result_plot,  
       height = 12, width = 24, units = "cm", dpi = 300)
