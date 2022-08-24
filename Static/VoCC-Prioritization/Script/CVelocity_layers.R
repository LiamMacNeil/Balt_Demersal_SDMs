library(tidync)
library(ncdf4)
library(ncmeta)
library(RNetCDF)
library(tidyverse)
library(lubridate)
library(timetk)
library(terra)
library(raster)
library(fields)
library(stars)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

physics_file <- "CMEMS-BALTICSEA_003_011-bottomT_mlotst_so_sob_thetao-2001_2020.nc"

physics <- tidync("CMEMS-BALTICSEA_003_011-bottomT_mlotst_so_sob_thetao-2001_2020.nc")

# extract datetime
time_ex_physics <- physics %>% activate("D0") %>% hyper_array()

tunit_phys <- ncmeta::nc_atts(physics_file, "time") %>% 
  #tidyr::unnest(cols = c(value)) %>% 
  dplyr::filter(name == "units")

time_parts_phys <- RNetCDF::utcal.nc(as.character(tunit_phys$value), time_ex_physics$time)

## convert to date-time
physics[["transforms"]][["time"]][["time"]] <- ISOdatetime(time_parts_phys[,"year"], 
                                                           time_parts_phys[,"month"], 
                                                           time_parts_phys[,"day"], 
                                                           time_parts_phys[,"hour"], 
                                                           time_parts_phys[,"minute"], 
                                                           time_parts_phys[,"second"])

# Surface variables are defined by median values 1-5m
physics_surface <- physics %>% activate(thetao) %>% 
  hyper_filter(depth = depth < 5) %>% 
  hyper_tibble() %>% 
  group_by(latitude, longitude, time) %>% 
  summarise_at(vars(thetao), funs(median(., na.rm = F)))

# format data in space-time
physics_surface <- physics_surface %>% arrange(time, longitude, latitude)

physics_surface$time <- year(physics_surface$time)

P_2001 <- physics_surface %>% 
  filter(time == 2001)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2001[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2001[, 4])

covs_vect <- vect(P_2001, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

r <- rast(extent = ext(covs_spatial),
          resolution = c(0.0575, 0.0575),
          vals=0,
          crs = as.character(st_crs(covs_spatial))[1])
temp <- raster(terra::rasterize(covs_vect, 
                                r,
                                "thetao"))
names(temp)<-"Surface temperature 2001"

### 2002 
P_2002 <- physics_surface %>% 
  filter(time == 2002)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2002[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2002[, 4])

covs_vect <- vect(P_2002, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

r <- rast(extent = ext(covs_spatial),
          resolution = c(0.0575, 0.0575),
          vals=0,
          crs = as.character(st_crs(covs_spatial))[1])
temp_2002 <- raster(terra::rasterize(covs_vect, 
                                r,
                                "thetao"))
names(temp_2002)<-"Surface temperature 2002"

### 2003 
P_2003 <- physics_surface %>% 
  filter(time == 2003)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2003[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2003[, 4])

covs_vect <- vect(P_2003, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2003 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2003)<-"Surface temperature 2003"

### 2004
P_2004 <- physics_surface %>% 
  filter(time == 2004)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2004[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2004[, 4])

covs_vect <- vect(P_2004, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2004 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2004)<-"Surface temperature 2004"

### 2005
###
###

P_2005 <- physics_surface %>% 
  filter(time == 2005)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2005[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2005[, 4])

covs_vect <- vect(P_2005, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2005 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2005)<-"Surface temperature 2005"

### 2006
###
###

P_2006 <- physics_surface %>% 
  filter(time == 2006)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2006[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2006[, 4])

covs_vect <- vect(P_2006, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2006 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2006)<-"Surface temperature 2006"

### 2007
###
###

P_2007 <- physics_surface %>% 
  filter(time == 2007)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2007[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2007[, 4])

covs_vect <- vect(P_2007, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2007 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2007)<-"Surface temperature 2007"

### 2008
###
###

P_2008 <- physics_surface %>% 
  filter(time == 2008)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2008[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2008[, 4])

covs_vect <- vect(P_2008, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2008 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2008)<-"Surface temperature 2008"

### 2009
###
###

P_2009 <- physics_surface %>% 
  filter(time == 2009)
# convert dataframe into spatial object for further analysis
covs_spatial <- SpatialPointsDataFrame(coords = P_2009[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2009[, 4])

covs_vect <- vect(P_2009, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")

temp_2009 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2009)<-"Surface temperature 2009"

### 2010
###
###

P_2010 <- physics_surface %>% 
  filter(time == 2010)
covs_spatial <- SpatialPointsDataFrame(coords = P_2010[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2010[, 4])
covs_vect <- vect(P_2010, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2010 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2010)<-"Surface temperature 2010"

### 2011
###
###

P_2011 <- physics_surface %>% 
  filter(time == 2011)
covs_spatial <- SpatialPointsDataFrame(coords = P_2011[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2011[, 4])
covs_vect <- vect(P_2011, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2011 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2011)<-"Surface temperature 2011"

### 2012
###
###

P_2012 <- physics_surface %>% 
  filter(time == 2012)
covs_spatial <- SpatialPointsDataFrame(coords = P_2012[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2012[, 4])
covs_vect <- vect(P_2012, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2012 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2012)<-"Surface temperature 2012"

### 2013
###
###

P_2013 <- physics_surface %>% 
  filter(time == 2013)
covs_spatial <- SpatialPointsDataFrame(coords = P_2013[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2013[, 4])
covs_vect <- vect(P_2013, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2013 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2013)<-"Surface temperature 2013"

### 2014
###
###

P_2014 <- physics_surface %>% 
  filter(time == 2014)
covs_spatial <- SpatialPointsDataFrame(coords = P_2014[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2014[, 4])
covs_vect <- vect(P_2014, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2014 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2014)<-"Surface temperature 2014"

### 2015
###
###

P_2015 <- physics_surface %>% 
  filter(time == 2015)
covs_spatial <- SpatialPointsDataFrame(coords = P_2015[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2015[, 4])
covs_vect <- vect(P_2015, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2015 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2015)<-"Surface temperature 2015"

### 2016
###
###

P_2016 <- physics_surface %>% 
  filter(time == 2016)
covs_spatial <- SpatialPointsDataFrame(coords = P_2016[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2016[, 4])
covs_vect <- vect(P_2016, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2016 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2016)<-"Surface temperature 2016"

### 2017
###
###

P_2017 <- physics_surface %>% 
  filter(time == 2017)
covs_spatial <- SpatialPointsDataFrame(coords = P_2017[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2017[, 4])
covs_vect <- vect(P_2017, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2017 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2017)<-"Surface temperature 2017"

### 2018
###
###

P_2018 <- physics_surface %>% 
  filter(time == 2018)
covs_spatial <- SpatialPointsDataFrame(coords = P_2018[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2018[, 4])
covs_vect <- vect(P_2018, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2018 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2018)<-"Surface temperature 2018"

### 2019
###
###

P_2019 <- physics_surface %>% 
  filter(time == 2019)
covs_spatial <- SpatialPointsDataFrame(coords = P_2019[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2019[, 4])
covs_vect <- vect(P_2019, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2019 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2019)<-"Surface temperature 2019"

### 2020
###
###

P_2020 <- physics_surface %>% 
  filter(time == 2020)
covs_spatial <- SpatialPointsDataFrame(coords = P_2020[, c("longitude", "latitude")], 
                                       proj4string = CRS("+proj=longlat +datum=WGS84"),
                                       data = P_2020[, 4])
covs_vect <- vect(P_2020, 
                  geom = c("longitude", "latitude"), 
                  crs = "+proj=longlat +datum=WGS84")
temp_2020 <- raster(terra::rasterize(covs_vect, 
                                     r,
                                     "thetao"))
names(temp_2020)<-"Surface temperature 2020"

temp_layers <- stack(temp,
                     temp_2002,
                     temp_2003, 
                     temp_2004,
                     temp_2005,
                     temp_2006,
                     temp_2007,
                     temp_2008, 
                     temp_2009,
                     temp_2010,
                     temp_2011,
                     temp_2012,
                     temp_2013,
                     temp_2014,
                     temp_2015,
                     temp_2016,
                     temp_2017,
                     temp_2018,
                     temp_2019,
                     temp_2020)

layers_dir <- paste0("SurfaceTemp_CMEMS")
dir.create(layers_dir)
for (l in 1:nlayers(temp_layers)) {
  writeRaster(temp_layers[[l]], filename = paste0(layers_dir, "/", names(temp_layers)[l], ".tif"))
}
