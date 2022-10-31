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

# Load all exchange data (unless already generated, then proceed below)

# Haul (HH) and catch data (CA)
HH_CA <- read_csv("../Data/ExchangeData_HH_CA_2022-09-21_15_30_03.csv")
HH_CA <- HH_CA %>% 
  filter(StNo != -9) %>% 
  mutate(Date=make_date(Year, Month, Day)) %>% 
  filter(Gear == "TVL"| Gear == "TVS") %>% 
  filter(HaulDur >! 90) %>% 
  filter(HaulDur >= 15) 
  
HH_CA_J <- HH_CA[,c("Quarter", "Gear", "StNo", 
               "HaulNo", "Year", "Date", "HaulDur",
               "ShootLat", "ShootLong", "HaulLat",
               "HaulLong", "Depth", "HaulVal", "Distance")]

# Lenth (HL)
HL <- read_csv("../Data/ExchangeData_HL_2022-09-21_15_30_03.csv")
HL <- HL %>% 
  filter(StNo != -9) %>% 
  filter(Gear == "TVL"| Gear == "TVS") %>% 
  filter(CatCatchWgt != "NA") %>% 
  filter(LngtClass != -9) 
    
HL <- HL[,c("Quarter", "Gear", "StNo", 
            "HaulNo", "Year", "Sex",
            "TotalNo", "CatCatchWgt", "LngtCode",
            "LngtClass", "ValidAphiaID", "ScientificName_WoRMS")]

# Date and StnNo matching

#Sweptarea_trim <- Sweptarea[,c("date", "StNo","HaulNo", "Year", "SweptAreaDSKM2", "meanLong", "meanLat")]

# Not applying "by" argument convieniently executes natural join by common columns
HH_HL_CA <- left_join(HL, HH_CA_J)

# Full data
HH_HL_CA <- HH_HL_CA %>% 
  filter(Year <= 2020)

#Cod and Herring
HH_HL_CA_Cod_Herring <- HH_HL_CA %>% 
  filter(ScientificName_WoRMS == "Gadus morhua" | ScientificName_WoRMS == "Clupea harengus") %>% 
  # Trim years to match Swept area
  filter(Year <= 2020)

Sweptarea <- read_csv("../Data/SweptAreaAssessmentOutput_2022-09-21 19_20_08.csv")
Sweptarea <- Sweptarea %>% 
  filter(Gear == "TVL" |Gear == "TVS") %>% 
  filter(SweptAreaDSKM2 != -9) %>% 
  mutate(date=make_date(Year, Month, Day))

HH_HL_CA_Cod_Herring_Swept <- left_join(HH_HL_CA_Cod_Herring, Sweptarea)
HH_HL_CA_AllSpec_Swept <- left_join(HH_HL_CA, Sweptarea)

HH_HL_CA_Cod_Herring_Swept <- HH_HL_CA_Cod_Herring_Swept %>% 
  mutate(Density_g = (CatCatchWgt/SweptAreaDSKM2)/1000) 

HH_HL_CA_AllSpec_Swept <- HH_HL_CA_AllSpec_Swept %>% 
  mutate(Density_g = (CatCatchWgt/SweptAreaDSKM2)/1000) 

write_csv(HH_HL_CA_Cod_Herring_Swept, "../Data/HH_HL_CA_Cod_Herring_Swept.csv")
write_csv(HH_HL_CA_AllSpec_Swept, "../Data/HH_HL_CA_AllSpec_Swept.csv")

# If already generated:

HH_HL_CA_Cod_Herring_Swept <- read_csv("../Data/HH_HL_CA_Cod_Herring_Swept.csv")
HH_HL_CA_AllSpec_Swept <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv")

# Plotting examples
Quarter <- c("Q1", "Q4")
names(Quarter) <- c(1, 4)

# Facet by species and color by quarter
HH_HL_CA_Cod_Herring_Swept %>%
  filter(Year != 2008) %>% 
  #filter(Year != 2002) %>% 
  ggplot(aes(Year, Density_g, colour = factor(Quarter))) +
  stat_summary(fun.data = "mean_cl_boot") +
  #expand_limits(y = 0) +
  scale_color_brewer(palette = "Set2") +
  guides(color=guide_legend(override.aes = list(size=1.5)))+
  labs(x = "Year", y = expression(Biomass~Density~(kg~km^2)),
       colour = "Quarter")+
  facet_wrap(~ ScientificName_WoRMS + Development)+
  theme_bw(24)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8))

# Multiple conditions column 
x <- HH_HL_CA_Cod_Herring_Swept %>% 
  #group_by(ScientificName_WoRMS) %>% 
  mutate(Development = ifelse((HH_HL_CA_Cod_Herring_Swept$ScientificName_WoRMS == "Gadus morhua" & LngtClass > 350 | HH_HL_CA_Cod_Herring_Swept$ScientificName_WoRMS == "Clupea harengus" & LngtClass > 250), 
                              "Adult", "Juvenile"))

# Facet grid by species x size and color by quarter
x %>% 
  filter(Year != 2008) %>% 
  ggplot(aes(Year, Density_g, colour = factor(Quarter))) +
  stat_summary(fun.data = "mean_cl_boot") +
  #expand_limits(y = 0) +
  scale_color_brewer(palette = "Set2") +
  guides(color=guide_legend(override.aes = list(size=1.5)))+
  labs(x = "Year", y = expression(Biomass~Density~(kg~km^2)),
       colour = "Quarter")+
  facet_grid( Development ~ ScientificName_WoRMS)+
  theme_bw(24)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8))

# Interpolation testing

balt <- rast("../../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
pr <- as.polygons(balt > -Inf)
pr <- as(pr, "Spatial")
plot(pr)
Baltic_sf <- st_as_sf(pr)

# map data
xlim <- c(9.25, 26.5)
ylim <- c(53.5, 58.5)

m <- map_data("worldHires", xlim = xlim, ylim = ylim)

HH_HL_CA_Cod_Herring_Swept <- read_csv("../Data/HH_HL_CA_Cod_Herring_Swept.csv")

# missing meanLat/meanLong values removed
Herring_Swept <- HH_HL_CA_Cod_Herring_Swept %>% 
  filter(meanLat != "NA") %>% 
  filter(meanLong != "NA") %>% 
  filter(ScientificName_WoRMS == "Clupea harengus") %>% 
  filter(Quarter == 4) %>% 
  filter(LngtClass >= 250) 

Cod_Swept <- HH_HL_CA_Cod_Herring_Swept %>% 
  filter(meanLat != "NA") %>% 
  filter(meanLong != "NA") %>% 
  filter(ScientificName_WoRMS == "Gadus morhua") %>% 
  filter(Quarter == 4) %>% 
  filter(LngtClass >= 350) 
  
dat_sf <- st_as_sf(Herring_Swept, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84")

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

# Using to estimate useful maxdist
plot(variogram(fit_idw))

grid$IDW <- predict(fit_idw, newdata = as(grid, "Spatial")) %>% 
  st_as_sf() %>% 
  pull(1)

raster <- grid %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(X,Y,IDW) %>% 
  rasterFromXYZ(crs=CRS("+proj=longlat +datum=WGS84"))

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

ggplot(data = HH_HL_CA_Cod_Herring_Swept, aes(meanLong, meanLat)) +
  geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
  coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
  geom_point(aes(color=Gear),size=0.3, alpha=0.1) +
  theme_bw(30)+
  theme(plot.title = element_text (hjust = 0.5)) +
  xlab(expression('Longitude '*degree*"W")) +
  ylab(expression("Latitude "*degree*"N")) +
  guides(color=guide_legend(override.aes = list(size=5, alpha = 1))) +
  scale_color_brewer(palette = "Dark2")+
  ggtitle("DATRAS-BITS Trawl 2001-2020")
  
  
# Fine tune power function through testing 
# https://mgimond.github.io/Spatial/interpolation-in-r.html#idw
library(stars)

# Leave-one-out validation routine
IDW.out <- vector(length = length(dat_sf))
for (i in 1:length(dat_sf)) {
  IDW.out[i] <- idw(Density_g ~ 1, dat_sf[-i,], dat_sf[i,], idp=2.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ dat_sf$Density_g, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ dat_sf$Density_g), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
