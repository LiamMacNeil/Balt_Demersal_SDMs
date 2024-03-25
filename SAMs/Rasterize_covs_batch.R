library(tidyverse)
library(lubridate)
library(raster)
library(stars)
library(sp)
library(terra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

covariates <- read_csv("../../../copernicus-processed-data/covariates/CMEMS_Layers_df/merged_covariates.csv")

icesarea_sf <- gisland::read_sf_ftp("ices_areas") %>%
  filter(subdivisio == "21" |subdivisio == "22" | subdivisio == "23" |subdivisio == "24"|subdivisio == "25"|subdivisio == "26"|subdivisio == "27"|subdivisio == "28"|subdivisio == "29") 

# Unit conversions for covariates
# Oxygen (mg/L)
covariates$o2 <- ((((covariates$o2)/1000)*31.998)/1000)*1000
covariates$o2b <- ((((covariates$o2b)/1000)*31.998)/1000)*1000

covariates <- covariates %>% 
  mutate(Year = year(time)) %>% 
  mutate(Quarter = quarter(time)) %>% 
  filter(Quarter == 1| Quarter == 4)
  

####################################
# Surface Temperature 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    temp_df <- filter(covariates, Year == i & Quarter == j) 
    temp_df_i <- temp_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, thetao) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(temp = median(thetao)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(temp_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "temp")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
  sfc_p <- filter(sfc, Year == i & Quarter == j)
  rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                             dplyr::select(geometry, temp))), "Raster")
  #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
  rast_resamp <- resample(rast, r, method = "bilinear")
  names(rast_resamp) <- paste0("thetao", "-",i, "-", j)
  stack <- raster::stack(stack, rast_resamp)
  }
  }

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                                         names(stack)[p], ".tif"))
}  

####################################
# Bottom Temperature 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    temp_df <- filter(covariates, Year == i & Quarter == j) 
    temp_df_i <- temp_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, bottomT) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(temp = median(bottomT)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(temp_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "temp")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, temp))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("bottomT", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  


####################################
# Surface Salinity 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    salin_df <- filter(covariates, Year == i & Quarter == j) 
    salin_df_i <- salin_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, so) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(salin = median(so)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(salin_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "salin")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, salin))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("so", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  

####################################
# Bottom Salinity 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    salin_df <- filter(covariates, Year == i & Quarter == j) 
    salin_df_i <- salin_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, sob) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(salin = median(sob)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(salin_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "salin")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, salin))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("sob", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  

####################################
# Surface Oxygen 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    o2_df <- filter(covariates, Year == i & Quarter == j) 
    o2_df_i <- o2_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, o2) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(oxy = median(o2)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(o2_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "oxy")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, oxy))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("o2", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  

####################################
# Bottom Oxygen 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    o2b_df <- filter(covariates, Year == i & Quarter == j) 
    o2b_df_i <- o2b_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, o2b) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(oxy = median(o2b)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(o2b_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "oxy")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, oxy))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("o2b", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  

####################################
# Chlorophyll 
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    chl_df <- filter(covariates, Year == i & Quarter == j) 
    chl_df_i <- chl_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, chl) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(chl = median(chl)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(chl_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "chl")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, chl))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("chl", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  

####################################
# MLD
####################################

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
final <- list()

# Run interpolation to points through sf
for(i in Year){
  for(j in Quarter){
    mlotst_df <- filter(covariates, Year == i & Quarter == j) 
    mlotst_df_i <- mlotst_df %>% 
      dplyr::select(longitude, latitude, time, Year, Quarter, mlotst) %>% 
      group_by(longitude, latitude, time, Year, Quarter) %>% 
      mutate(mlotst = median(mlotst)) %>% 
      ungroup() %>% 
      as_tibble()
    #temp_sf_i <- st_as_sf(temp_df_i, coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84")
    final <- append(final, list(mlotst_df_i[,c("longitude", "latitude", "time", "Year", "Quarter", "mlotst")]))
  }
}

# bind into sf lists
sfc = do.call(bind_rows, final)

Year <- unique(covariates$Year)
Quarter <- unique(covariates$Quarter)
stack <- raster::stack()

r <- raster("../Predictors/CMEMS_Habitat_cut/Surface.temperature.tif")
r <- crop(mask(r, icesarea_sf), extent(icesarea_sf))

for(i in Year){
  for(j in Quarter){  
    sfc_p <- filter(sfc, Year == i & Quarter == j)
    rast <- as(st_rasterize((sfc_p %>% st_as_sf(coords = c("longitude", "latitude")) %>% 
                               dplyr::select(geometry, mlotst))), "Raster")
    #rast_crop <- crop(mask(rast, icesarea_sf), extent(icesarea_sf))
    rast_resamp <- resample(rast, r, method = "bilinear")
    names(rast_resamp) <- paste0("mlotst", "-",i, "-", j)
    stack <- raster::stack(stack, rast_resamp)
  }
}

for (p in 1:nlayers(stack)) {
  writeRaster(stack[[p]], 
              filename = paste0("../Predictors/Q1-Q4_full/",
                                names(stack)[p], ".tif"))
}  