library(tidyverse)
library(sf)
library(oce)
library(raster)
library(stars)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

xlim <- c(9.420922, 22.51811)
ylim <- c(53.91644, 59.23877)

#To be deal some stubbornly invalid geometries in the shapefiles 
sf_use_s2(F)

# High res polygons
icesarea <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

Coastline <- read_sf("../../Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
  st_transform(4326) %>% 
  st_crop(st_bbox(icesarea))

taxa <- c("Pleuronectes platessa", "Platichthys flesus", 
          "Adult Gadus morhua","Juvenile Gadus morhua","Limanda limanda")

# Raw quantities
dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv") %>% 
  filter(HaulVal == "V" | HaulVal == "C"| HaulVal == "A") %>% 
  filter(meanLong != "NA") %>% 
  filter(meanLat != "NA") %>% 
  mutate(Depth = abs(Depth)) %>% 
  mutate(Density_kg = Density_g/1000 ,
         Density_kg_log10 = log(Density_kg+1)) %>% 
  mutate(Haul_ID = paste0(Survey, "_", Year, "_", Quarter, 
                          "_", Country, "_", Ship, "_", Gear, 
                          "_", Depth)) %>% 
  dplyr::select(Quarter, Gear, HaulNo, Year, CatCatchWgt,
                Depth, meanLat, meanLong, Flag_SweptAreaDSKM2, date,
                Density_kg,Density_kg_log10,ScientificName_WoRMS,LngtClass,Haul_ID) %>% 
  mutate(Species_group = factor(case_when(LngtClass > 34 & ScientificName_WoRMS == 'Gadus morhua' ~ "Adult Gadus morhua",
                                          LngtClass < 35 & ScientificName_WoRMS == 'Gadus morhua' ~ "Juvenile Gadus morhua",
                                          TRUE ~ ScientificName_WoRMS))) %>% 
  group_by(Haul_ID, Species_group, date, Year,Gear,
           Depth, meanLat, meanLong, Quarter, .drop = FALSE) %>% 
  summarise(Density_kg = sum(Density_kg)) %>% 
  ungroup() %>% 
  fill(Haul_ID, Species_group, date, Year,
       Depth, Quarter, meanLat, meanLong,Gear) %>% 
  arrange(date, meanLat, meanLong) %>% 
  filter(Species_group %in% taxa)  %>% 
  st_as_sf(coords = c("meanLong", "meanLat"), 
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea)

# Extra for cropping tiled sf smooths
###########################
icesarea_sf <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "20" |SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

BITS_buffer <- dat %>% 
  #st_as_sf(coords = c("Longitude", "Latitude"), 
  #         crs = st_crs(4326)) %>% 
  st_buffer(dist = 0.5) %>% 
  st_geometry() %>% 
  st_union() %>% 
  st_crop(xmin=9.420922, ymin=53.91644, xmax=22.51811, ymax=59.23877)

############################
#
############################
gebco <- raster("../../Oceanographic/Data/gebco_2022_sub_ice_topo_geotiff/gebco_2022_sub_ice_n90.0_s0.0_w0.0_e90.0.tif")
highres_bathy_rast <- crop(mask(gebco, icesarea), 
                           extent(icesarea))
#bathy_stars <- st_as_sf(st_as_stars(highres_bathy_rast))
#deeps <- highres_bathy_rast < -70
#deeps_sf <- st_as_sf(st_as_stars(deeps)) %>% 
#  filter(layer == TRUE) %>% 
#  st_union() 
rast_spdf <- as.data.frame(as(highres_bathy_rast, "SpatialPixelsDataFrame"))
colnames(rast_spdf) <- c("depth", "x", "y")
rast_spdf <- rast_spdf %>% 
  filter(depth < 0.1) %>% 
  mutate(depth = depth*-1)

##########################
# Baltic
##########################
deeps <- colormap(colormaps$velocity_blue, nshades = 10)

balt <- ggplot() +
  geom_tile(data=rast_spdf, aes(x=x, y=y, fill=depth)) + 
  geom_sf(data=dat, aes(color=Gear),  size=0.5, alpha=0.01)+
  #geom_sf(data = BITS_buffer, fill=NA, linewidth=0.6, linetype=2, col = "black")+
  geom_sf(data = icesarea_sf, fill=NA, linewidth=0.25, linetype=1, col = "black")+
  geom_sf(data = Coastline, fill="grey70", color="black", linewidth=0.05)+
  theme_classic(16)+
  coord_sf(expand = F)+
  xlab(expression('Longitude'))+
  ylab(expression("Latitude"))+
  scale_fill_steps(#low = deeps[10],
                   #high = deeps[1],
                   low = oce.colorsGebco(10)[10],
                   high = oce.colorsGebco(10)[1],
                   name = "Depth (m)",
                   n.breaks = 10,
                   breaks = seq(0,400, 50))+
  scale_color_manual(values = c("#009E73", "#882255"),
                     name = "Gear Type", 
                     labels = c("TVL (930 meshes)",
                                "TVS (520 meshes)")) +
  guides(fill = guide_legend(override.aes = list(size=2)),
         color = guide_legend(override.aes = list(size = 5, alpha = 1))) 
ggsave("../Figures/Spring2024Revision/BITS_TrawlCoverage.png", balt, 
       width = 17, height = 15, units = "cm", dpi = 600)

##########################
# Europe
##########################
icesarea_sf <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29"|
           SubDivisio == "30"|SubDivisio == "31"|SubDivisio == "32") %>% 
  st_transform(st_crs(4326)) 

### High res polygons
LMEs_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/LMEs/", layer = "LMEs66") %>% 
  filter(LME_NAME == "Baltic Sea" | LME_NAME == "North Sea" | LME_NAME == "Celtic-Biscay Shelf" | 
           LME_NAME == "Iberian Coastal" | LME_NAME == "Norwegian Sea" |
           LME_NAME == "Barents Sea" ) %>% 
  mutate(Coast = case_when(LME_NAME == "Baltic Sea" ~ "Europe",
                           LME_NAME == "North Sea" ~ "Europe",
                           LME_NAME == "Celtic-Biscay Shelf" ~ "Europe",
                           LME_NAME == "Iberian Coastal" ~ "Europe",
                           LME_NAME == "Norwegian Sea" ~ "Europe",
                           LME_NAME == "Barents Sea" ~ "Europe")) %>% 
  st_transform(3995) 

LME_names <- as.data.frame(LMEs_shelf_seas[,c("LME_NUMBER",
                                              "LME_NAME",
                                              "Coast")])[,c(1,2,3,4)]

Ecoregions_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/MEOW/", layer = "meow_ecos")%>% 
  st_transform(3995) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))

Coastline <- read_sf("../../Oceanographic/Data/GSHHS_shp/l/", layer = "GSHHS_l_L1") %>% 
  st_transform(3995) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))

europe <- ggplot()+
  geom_sf(data=Ecoregions_shelf_seas, fill=NA, color=NA, linewidth=1e-5)+
  #geom_sf(data = BITS_buffer, fill=NA, color="black",linewidth=0.5, linetype=2)+
  geom_sf(data=Coastline,fill="grey70", color="black", linewidth=0.05)+
  geom_sf(data = icesarea_sf, fill=NA, linewidth=0.1, linetype=1, col = "#D55E00")+
  theme_classic(16)+
  coord_sf(expand = F)+
  labs(x="", y="")+
  theme(axis.text = element_blank())
  
ggsave("../Figures/Spring2024Revision/Europe_TrawlCoverage.png", europe, 
       width = 17, height = 15, units = "cm", dpi = 600)
