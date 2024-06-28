library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(raster)
library(viridis)
library(itsadug)
library(oce)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#data
dat <- read_csv("../Data/Taxa_env_GAMs_v2_cutcovs.csv") %>% 
  rename(SurfaceOxygen = o2...5, 
         BottomOxygen = o2...6,
         SurfaceSalinity = so...7,
         BottomSalinity = so...8
  ) %>% 
  mutate(ScientificName_WoRMS = factor(ScientificName_WoRMS, 
                                       levels = c("Limanda limanda",
                                                  "Platichthys flesus",
                                                  "Pleuronectes platessa",
                                                  "Juvenile Gadus morhua",
                                                  "Adult Gadus morhua")),
         Quarter = factor(Quarter),
         Density_log = log(Density_kg+1)) %>% 
  mutate(Depth = abs(Depth)) %>% 
  filter(Depth < 121) %>% 
  #filter(Year !=2008 | Density_log < 7) %>% 
  #filter(Year != 2008 & 
  #         ScientificName_WoRMS != "Juvenile Gadus morhua" | 
  #         ScientificName_WoRMS != "Adult Gadus morhua") %>% 
  mutate(Year_fac = factor(Year, levels =c(2001:2020), ordered = T)) %>% 
  as.data.frame()

xlim <- c(9.420922, 22.51811)
ylim <- c(53.91644, 59.23877)

#To be deal some stubbornly invalid geometries in the shapefiles 
sf_use_s2(F)

# High res polygons
icesarea <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) %>% 
  as("Spatial")  

icesarea_sf <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

Coastline <- read_sf("../../Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
  st_transform(4326) %>% 
  st_crop(st_bbox(icesarea))

#icesarea_sf <- read_sf("../../../Oceanographic/Data/ICES_shapefiles/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
#  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
#           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
#  st_transform(st_crs(4326)) %>% 
#  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

# Extra for cropping tiled sf smooths
###########################
BITS_buffer <- dat %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = st_crs(4326)) %>% 
  st_buffer(dist = 0.5) %>% 
  st_geometry() %>% 
  st_union()

# Quarterly, Decadal Predictors
Layers <- stack(list.files("../Predictors/Yearly_quarterly_predcitors", 
                           pattern = "\\.tif$", full.names = TRUE))
Layers_cut <- crop(mask(Layers,as(BITS_buffer, "Spatial")), extent(as(BITS_buffer, "Spatial")))[[1]]


GAM_3 <- readRDS("../GAMs/GAM3.rds")

# how to spatially smooth and completely predict?

dat$pred <- predict(GAM_3, newdata = dat %>% 
                      dplyr::select(Latitude, Longitude, Quarter,
                                    Year_fac, ScientificName_WoRMS, 
                                    Depth,Density_log, bottomT,
                                    BottomOxygen, BottomSalinity), 
                    type="response")

#################################################################
# Spatially gridding
#################################################################

sf_points <- st_as_sf(dat, 
                      coords = c("Longitude", "Latitude"), 
                      crs = 4326) 
index <- which(lengths(st_intersects(BITS_buffer, sf_points))>0)

fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}
#############################################
# Summed grid

full_grid <- fishnet(sf_points, square = F, cellsize = c(0.4, 0.4))

full_output <- full_grid %>% 
  st_as_sf(wkt = "geometry") %>% 
  mutate(grid_id = row_number()) %>% 
  st_join(sf_points, left=T) %>% 
  arrange(Year)

#############################################


#Years <- unique(sf_points$Year)
#dynamic <-list()

#for(i in Years){
  
 # x <- sf_points %>% filter(Year == i)
  #grid <- fishnet(x, square = F, cellsize = c(0.5, 0.5))
  
  #output <- grid %>% 
  #  st_as_sf(wkt = "geometry") %>% 
  #  mutate(grid_id = row_number()) %>% 
  #  st_join(x)
  
  #dynamic[[i]] <- output
#}

#Allyears <- do.call("rbind", dynamic)

#Allyears_gridID <- Allyears %>% 
#  rename(geometry = x.x) %>% 
#  complete(grid_id, pred) %>% 
#  st_as_sf() 
  
map <- full_output %>% 
  filter(ScientificName_WoRMS == "Limanda limanda") %>% 
  ggplot()+
  geom_sf(aes(fill = pred), color = "black", linewidth=0.01)+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf(expand = F)+
  scale_fill_gradientn(colors = oce.colorsDensity(20),
                       name = expression(Predicted~log(kg~km^-2)),
                       na.value = NA)+
  theme_classic(10)+
  theme(axis.text.x.bottom = element_text(angle=30, hjust=1),
        axis.text = element_text(size=6))
ggsave("../Figures/Spring2024Revision/GAM3_Dab_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

#################################################################
# for each species, for each grid cell, calculate early (2001-2005)-later (2016-2020) averages
early <- c(2001:2005)
late <- c(2016:2020) 

delta <- full_output %>%
  rename(geometry = x.x) %>% 
  mutate(Phase = case_when(Year %in% early ~ "Before",
                           Year %in% late ~ "After",
                           .default = "FALSE")) %>% 
  filter(Phase != "FALSE") %>% 
  group_by(ScientificName_WoRMS, grid_id, Phase) %>% 
  summarise(mean = mean(pred)) %>% 
  ungroup() %>% 
  complete(ScientificName_WoRMS, grid_id, Phase) %>% 
  pivot_wider(names_from = Phase, values_from = mean) %>% 
  drop_na() %>% 
  mutate(diff = After - Before) %>% 
  st_as_sf()


#coords <-   full_output %>% 
#  rename(geometry = x.x) %>% 
#  distinct(grid_id, geometry)
#delta <- delta %>% 
  #full_join(coords) %>% 
  #mutate(diff = f-l)

limit <- max(abs(delta$diff)) * c(-1, 1)

map <- delta %>%
  mutate(ScientificName_WoRMS = factor(ScientificName_WoRMS, 
                                       levels = c("Limanda limanda",
                                                  "Platichthys flesus",
                                                  "Pleuronectes platessa",
                                                  "Juvenile Gadus morhua",
                                                  "Adult Gadus morhua"),
                                       labels = c("Dab",
                                                  "Flounder",
                                                  "Plaice",
                                                  "Juvenile Cod",
                                                  "Adult Cod"))) %>% 
  ggplot()+
  geom_sf(aes(fill = diff), color = "black", linewidth=0.05)+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.25)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~ScientificName_WoRMS, nrow=3)+
  coord_sf(expand = F)+
  scale_fill_distiller(palette = "RdBu",
                       name = expression(log(kg~km^-2)),
                       na.value = NA, 
                       limit=limit,
                       direction = 1)+
  theme_bw(14)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1),
        axis.text = element_text(size=9))
ggsave("../Figures/Spring2024Revision/GAM3_pred_diffs.png", map,
       width = 17, height = 17, units = "cm", dpi = 600)
  

