library(raster)
library(terra)
library(mapdata)
library(tidyverse)
library(gisland)
library(sf)
library(gstat)
library(oce)
library(lubridate)
library(automap)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv")
dat$TotalNo <- as.numeric(dat$TotalNo)

mains = c("Gadus morhua", "Clupea harengus", "Pleuronectes platessa",
          "Platichthys flesus", "Scophthalmus maximus")

# missing meanLat/meanLong values removed
dat <- dat %>% 
  filter(meanLat != "NA") %>% 
  filter(meanLong != "NA") %>% 
  #filter(ScientificName_WoRMS == "Gadus morhua") %>% 
  #filter(ScientificName_WoRMS == "Clupea harengus") %>% 
  #filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  #filter(ScientificName_WoRMS == "Pleuronectes platessa") %>% 
  #filter(ScientificName_WoRMS == "Scophthalmus maximus") %>% 
  # For panel plotting
  filter(ScientificName_WoRMS %in% mains)  
  #filter(Depth > 20) 
  #mutate(Indkm2 = TotalNo/SweptAreaDSKM2)

dat %>% 
  filter(Year != 2008) %>% 
  ggplot(aes(Year, Density_g, colour = factor(Quarter))) +
  stat_summary(fun.data = "mean_cl_boot") +
  #expand_limits(y = 0) +
  scale_color_brewer(palette = "Set2") +
  guides(color=guide_legend(override.aes = list(size=1.5)))+
  labs(x = "Year", y = expression(Biomass~Density~(kg~km^2)),
       colour = "Quarter")+
  facet_grid(~ScientificName_WoRMS, scales = "free_y")+
  theme_bw(24)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8))

# Auxillary data
# map data
xlim <- c(9.25, 26.5)
ylim <- c(53.5, 58.5)
m <- map_data("worldHires", xlim = xlim, ylim = ylim)

balt <- rast("../../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
pr <- as.polygons(balt > -Inf)
pr <- as(pr, "Spatial")
plot(pr)
Baltic_sf <- st_as_sf(pr)

# Simple feature (sf) arrangement
dat_sf <- st_as_sf(dat, coords = c("meanLong", "meanLat"), crs = "+proj=longlat +datum=WGS84")

ggplot() + 
  geom_sf(data = Baltic_sf, fill = "grey") +
  geom_sf(data = dat_sf ,aes(fill = (Indkm2)), size=0.001)+
  theme_bw(12)+
  facet_grid(ScientificName_WoRMS ~ Quarter)+
  xlab(expression('Longitude '*degree*"E"))+
  ylab(expression("Latitude "*degree*"N"))+
  scale_fill_gradientn(colours = oce.colorsViridis(10), 
                       name=expression(kg~km^2))+  
  coord_sf()


# Interpolation testing
boundary_points <- st_cast(dat_sf, "POINT") %>% 
  #Wgt_sf <- rbind(boundary_points, Wgt_sf) %>%
  cbind(., st_coordinates(.))

grid <- st_make_grid(dat_sf, cellsize = c(0.2, 0.2), what = "centers") %>% 
  st_as_sf() %>% 
  cbind(.,st_coordinates(.))

fit_idw <- gstat(
  formula = SweptAreaDSKM2 ~ 1,
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
  theme_bw(30) +
  xlab(expression('Longitude '*degree*"E"))+
  ylab(expression("Latitude "*degree*"N"))+
  scale_fill_gradientn(colours = oce.colorsTemperature(10), 
                       #name=expression(kg~km^2))
                       #name=expression(Ind~km^2))
                       name=expression(km^2), n.breaks = 5)
                       



result = list()

for(i in mains) {
  f = as.formula(Density_g ~ 1)
  dat_sp <- filter(dat_sf, ScientificName_WoRMS == i)
  #v = autofitVariogram(f, as(dat_sf, "Spatial"))
  g = gstat(formula = f, data = as(dat_sp, "Spatial"), maxdist = 50, set = list(idp = 0.5))
  z = predict(g, as(grid, "Spatial"))
  names(z) = i
  z = z[,i]
  result[[i]] = z
}

result$along = 3
result = do.call(c, result)

a <- st_as_sf(unlist(result$`Gadus morhua`))
b <- st_as_sf(unlist(result$`Clupea harengus`))
c <- st_as_sf(unlist(result$`Pleuronectes platessa`))
d <- st_as_sf(unlist(result$`Platichthys flesus`))
e <- st_as_sf(unlist(result$`Scophthalmus maximus`))
