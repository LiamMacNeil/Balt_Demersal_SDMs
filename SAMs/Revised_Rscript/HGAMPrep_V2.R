library(tidyverse)
library(terra)
library(raster)
library(dismo)
library(sf)
library(mgcv)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################## ###################### ####################### #######################     
# Loading Data
####################### ####################### ####################### #######################
xlim <- c(9.420922, 22.51811)
ylim <- c(53.91644, 59.23877)

sf_use_s2(F)
icesarea <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) %>% 
  as("Spatial")  

# Quarterly, Decadal Predictors
Layers <- stack(list.files("../Predictors/Yearly_quarterly_predcitors", 
                           pattern = "\\.tif$", full.names = TRUE))
Layers_cut <- crop(mask(Layers,icesarea), extent(icesarea))

Coastline <- read_sf("../../../Feb2023_Transfer/Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
  st_transform(4326) %>% 
  st_crop(st_bbox(icesarea))

taxa <- c("Pleuronectes platessa", "Platichthys flesus", 
          "Adult Gadus morhua","Juvenile Gadus morhua","Limanda limanda")

# Raw quantities
dat <- read_csv("../Data/HH_HL_CA_AllSpec_Swept.csv") %>% 
  filter(HaulVal == "V" | HaulVal == "C"| HaulVal == "A") %>% 
  #filter(ScientificName_WoRMS %in% taxa) %>% 
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
                                          TRUE ~ ScientificName_WoRMS))) 

# Create trawl level data to model and map trawl level indices
dat_trawl <- dat %>% 
  group_by(Haul_ID, Species_group, date, Year,
           Depth, meanLat, meanLong, Quarter, .drop = FALSE) %>% 
  summarise(Density_kg = sum(Density_kg)) %>% 
  ungroup() %>% 
  fill(Haul_ID, Species_group, date, Year,
       Depth, Quarter, meanLat, meanLong) %>% 
  arrange(date, meanLat, meanLong) %>% 
  filter(Species_group %in% taxa)  
  

###########################################################
# Extracting environmental data
###########################################################

ind <- list()
Years <- 2001:2020
Quarters <- c(1,4)

Groups <- unique(dat_trawl$Species_group)
#Layers_cut[[which(grepl(paste0(i, ".", j), names(Layers_cut))==TRUE)]]

for(i in Years){
  #for(j in Quarters){
    
    t <- filter(dat_trawl, Year == i & Quarter == Quarters[2] & Species_group == Groups[5])
    covs <- Layers_cut[[which(grepl(paste0(i, ".", Quarters[2]), names(Layers_cut))==TRUE)]]
    
    spatDat <- SpatialPointsDataFrame(t[,c("meanLong","meanLat")], 
                                    data = data.frame(t[,c("Density_kg",
                                                           "Species_group",
                                                           "Depth")]))
    
    rast <- raster::rasterize(spatDat, 
                            covs[[1]], 
                            field="Density_kg")
    pts <- rasterToPoints(rast,
                        #,
                        #fun = function(x){
                        #  x > 0
                        #}
                        )
    rast_depth <- raster::rasterize(spatDat, 
                                    covs[[1]], 
                                    field="Depth")

    depth_pts <- rasterToPoints(rast_depth
                              #,
                              #fun = function(x){
                              #  x > 0
                              #}
                              )
    presence_cov <- raster::extract(covs, pts[,1:2])
    presence_cov <- data.frame(presence_cov)
    depth_comb <- data.frame(depth_pts)
    
    presence_cov$Density_kg <- pts[,3]
    presence_cov$Longitude <- pts[,1]
    presence_cov$Latitude <- pts[,2]
    
    depth_comb$Longitude <- depth_pts[,1]
    depth_comb$Latitude <- depth_pts[,2]
    
    full <- presence_cov
    
    full <- full %>% 
      full_join(depth_comb) %>% 
      rename(Depth = layer)
    
    full <- full[complete.cases(full),]
    full <- full[!is.infinite(rowSums(full)),]
    
    full$Year <- i
    full$Quarter <- Quarters[2]
    full$ScientificName_WoRMS <- Groups[5]
    
    colnames(full)[grepl('bottomT',colnames(full))] <- 'bottomT'
    colnames(full)[grepl('thetao',colnames(full))] <- 'temp'
    colnames(full)[grepl('chl',colnames(full))] <- 'chl'
    colnames(full)[grepl('mlotst',colnames(full))] <- 'mlotst'
    colnames(full)[grepl('o2',colnames(full))] <- 'o2'
    colnames(full)[grepl('o2.1',colnames(full))] <- 'o2b'
    colnames(full)[grepl('so',colnames(full))] <- 'so'
    colnames(full)[grepl('so.1',colnames(full))] <- 'sob'
    
    ind[[i]] <- full
#}
}

AdultCod_Q1 <- do.call("rbind", ind)
AdultCod_Q4 <- do.call("rbind", ind)

JuvenileCod_Q1 <- do.call("rbind", ind)
JuvenileCod_Q4 <- do.call("rbind", ind)

Dab_Q1 <- do.call("rbind", ind)
Dab_Q4 <- do.call("rbind", ind)

Flounder_Q1 <- do.call("rbind", ind)
Flounder_Q4 <- do.call("rbind", ind)

Plaice_Q1 <- do.call("rbind", ind)
Plaice_Q4 <- do.call("rbind", ind)



st <- do.call("rbind", list(Plaice_Q1,Plaice_Q4,
                            JuvenileCod_Q1,JuvenileCod_Q4,
                            AdultCod_Q1, AdultCod_Q4,
                            Flounder_Q1,Flounder_Q4, 
                            Dab_Q1,Dab_Q4))

write.csv(st, "../Data/Taxa_env_GAMs_v2.csv")

