library(icesDatras)
library(icesVocab)
library(raster)
library(terra)
library(mapdata)
library(tidyverse)
library(gisland)
library(sf)
library(gstat)
library(oce)
library(lubridate)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

balt <- rast("../../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
pr <- as.polygons(balt > -Inf)
pr <- as(pr, "Spatial")
plot(pr)
Baltic_sf <- st_as_sf(pr)

# map data
xlim <- c(9.25, 26.5)
ylim <- c(53.5, 58.5)

m <- map_data("worldHires", xlim = xlim, ylim = ylim)

icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  as("Spatial")

spatDat <- SpatialPointsDataFrame(Sweptarea[,c("meanLong","meanLat")], data.frame(Sweptarea[,"SweptAreaWSKM2"]))
min <- raster::rasterize(spatDat, layers_cut[[1]], field="SweptAreaWSKM2")

Sweptarea <- read_csv("../Data/SweptAreaAssessmentOutput_2022-09-21 19_20_08.csv")
  
Sweptarea <- Sweptarea %>% 
  filter(Gear == "TVL" |Gear == "TVS") %>% 
  filter(SweptAreaWSKM2 != "NA")
  
Swept_Q1 <- Sweptarea %>% 
  filter(Quarter == 1)
Swept_Q1 <- Sweptarea %>% 
  filter(Quarter == 4)

Swept_Q1_sf <- st_as_sf(Swept_Q1, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84") 
Swept_Q4_sf <- st_as_sf(Swept_Q1, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84") 

Swept_Q1_sf %>% ggplot() + 
  geom_sf(data = Baltic_sf, fill="grey80") +
  geom_sf(data = Swept_Q1_sf, aes(color = SweptAreaDSKM2), size=0.25) +
  theme_bw(20)+
  coord_sf()+
  ggtitle("Q4 Swept Area Coverage")+
  guides(color = guide_legend(override.aes = aes(size=4)))

SweptFig  <- Sweptarea %>% 
  ggplot(aes(x = Year, y = SweptAreaWSKM2)) +
  geom_col()+
  theme_bw(20)+
  labs(y = expression(Swept~Area~(km^2)))

ggsave("../Figures/SweptArea_2001_2020.png",SweptFig, height = 14, width = 17, units = "cm", dpi = 300)

########### 
########### Interpolated Versions
########### 

dat <- read_csv("../Data/Cod_Herring_CatchWgt_Density.csv")
dat <- dat %>% 
  filter(Taxa == "Herring") %>% 
  mutate(Density_g = Density/1000)
#%>% 
  #filter(Year.x !=2008)

dat_sf <- st_as_sf(dat, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84") 

dat <- read_csv("../Data/Cod_Herring_CatchWgt_Density_doorswept.csv")
dat <- dat %>% 
  filter(Taxa == "Cod") %>% 
  mutate(Density_g = Density/1000)

dat_sf <- st_as_sf(dat, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84")

boundary_points_Q1 <- st_cast(Swept_Q1_sf, "POINT") %>% 
  #Wgt_sf <- rbind(boundary_points, Wgt_sf) %>%
  cbind(., st_coordinates(.))
boundary_points_Q4 <- st_cast(Swept_Q4_sf, "POINT") %>% 
  #Wgt_sf <- rbind(boundary_points, Wgt_sf) %>%
  cbind(., st_coordinates(.))

boundary_points <- st_cast(dat_sf, "POINT") %>% 
  #Wgt_sf <- rbind(boundary_points, Wgt_sf) %>%
  cbind(., st_coordinates(.))

grid <- st_make_grid(dat_sf, cellsize = c(0.2, 0.2), what = "centers") %>% 
  st_as_sf() %>% 
  cbind(.,st_coordinates(.))

fit_idw <- gstat(
  formula = Density_g ~ 1,
  data = as(dat_sf, "Spatial"),
  maxdist = 50,
  set = list(idp = 0.5)
)    

################# Variogram 
plot(variogram(fit_idw))

grid$IDW <- predict(fit_idw, newdata = as(grid, "Spatial")) %>% 
  st_as_sf() %>% 
  pull(1)

raster <- grid %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(X,Y,IDW) %>% 
  rasterFromXYZ(crs=CRS("+proj=longlat +datum=WGS84"))

#raster<- crop(mask(raster, icesarea), extent(icesarea))

coords <- xyFromCell(raster[[1]], seq_len(ncell(raster[[1]])))
df <- stack(as.data.frame(getValues(raster[[1]])))
names(df) <- c('value', 'variable')
df$variable <- "Predicted"
df <- cbind(coords, df)

df %>% 
  ggplot() +
  geom_tile(aes(x, y, fill = value))+
  geom_polygon(data = m, aes(long, lat, group = group), fill = "white") +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  #geom_tile(aes(x, y, fill = value), data = df, inherit.aes = T)+
  #geom_path(data = icesarea, aes(long, lat, group = group), size=1) +
  theme_bw(18) +
  xlab(expression('Longitude '*degree*"E"))+
  ylab(expression("Latitude "*degree*"N"))+
  scale_fill_gradientn(colours = oce.colorsViridis(10), name=expression(kg~km^2))
  #scale_fill_gradientn(colours = oce.colorsViridis(10), name=expression(Swept~Area~(km^2)))

ggsave("../Figures/Q4_SweptArea_interp.svg", Swept_area_fig, height = 22, width = 16, units = "cm",
       dpi = 300)

### Catch Weight Modelling
aphia <- icesVocab::findAphia(c("Gadus morhua", "Clupea harengus"), latin = T)

Totalwgt <- getCatchWgt(survey = "BITS", years = 2001:2020, quarters = c(1,4), aphia = aphia)

Totalwgt$Taxa <- Totalwgt$Valid_Aphia
Totalwgt$Taxa <- ifelse(Totalwgt$Taxa == 126436, "Cod", "Herring") 

dat <- dat %>% 
  filter(Taxa == "Herring") %>% 
  mutate(Density_g = Density/1000)

dat <- 
  Totalwgt %>% 
  #filter(Distance != "NA") %>% 
  filter(HaulDur >! 90) %>% 
  filter(HaulDur >= 15) %>% 
  filter(Gear == "TVL"| Gear == "TVS") %>% 
  filter(CatchWgt != "NA") %>% 
  filter(HaulVal == "V" | HaulVal == "M")

# Create matching intercompatible dates
Sweptarea <- Sweptarea %>% 
  mutate(date=make_date(Year, Month, Day))
dat <- dat %>% 
  mutate(date=make_date(Year, Month, Day))
Sweptarea_trim <- Sweptarea[,c("date", "StNo","HaulNo", "Year", "SweptAreaDSKM2", "meanLong", "meanLat")]

dat_wgt_swept <- left_join(Sweptarea_trim, dat, by=c("StNo", "date"))

dat_wgt_swept <- dat_wgt_swept %>% 
  mutate(Density = CatchWgt / SweptAreaDSKM2) %>% 
  filter(Density != "NA")

Quarter <- c("Q1", "Q4")
names(Quarter) <- c(1, 4)

p <- dat_wgt_swept %>% 
  filter(Quarter !="NA") %>% 
  filter(Taxa == "Herring") %>% 
  ggplot() +
  geom_polygon(data = m, aes(long, lat, group = group), fill = "grey50") +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  geom_point(aes(x=meanLong, y=meanLat,
                 color=Gear), size=0.001)+
  facet_wrap(~Quarter,
             labeller = labeller(Quarter=Quarter))+
  theme_bw(24)+
  guides(color=guide_legend(override.aes = list(size=3)))+
  xlab(expression('Longitude '*degree*"E"))+
  ylab(expression("Latitude "*degree*"N"))+
  ggtitle("BITS Trawl Location 2001-2021")  

ggsave("../Figures/Quarterly_Trawl_Coverage_QAQC.svg", width = 22, height = 20, units = "cm",
       dpi = 300)

# 2008 Has enormous density values, remove them (for now)
dat_wgt_swept %>% filter(Year.x != 2008) %>% 
  ggplot(aes(x = Year.x, y = Density)) +
  geom_point()

p1 <- dat %>% 
  filter(Year.x != 2008) %>% 
  filter(Year.x != 2002) %>% 
  mutate(Density_g = Density/1000) %>% 
  ggplot(aes(Year.x, Density_g, colour = factor(Quarter))) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.6) +
  expand_limits(y = 0) +
  scale_color_brewer(palette = "Set2") +
  guides(color=guide_legend(override.aes = list(size=1.5)))+
  labs(x = "Year", y = expression(Biomass~Density~(kg~km^2)),
       colour = "Quarter")+
  facet_wrap(~ Taxa, scale = "free_y")+
  theme_bw(24)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8))
ggsave("../Figures/BiomassDensity_CodHerring.svg", p1, dpi = 300, height = 24, width = 35, units = "cm")

write_csv(dat_wgt_swept,"../Data/Cod_Herring_CatchWgt_Density_doorswept.csv" )


