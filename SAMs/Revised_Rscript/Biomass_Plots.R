library(tidyverse)
library(sf)

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

Coastline <- read_sf("../../../Feb2023_Transfer/Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
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
         Density_kg_log10 = log(Density_kg+1),
         CPUE = (as.numeric(TotalNo)/HaulDur)*60) %>% 
  mutate(Haul_ID = paste0(Survey, "_", Year, "_", Quarter, 
                          "_", Country, "_", Ship, "_", Gear, 
                          "_", Depth)) %>% 
  dplyr::select(Quarter, Gear, HaulNo, Year, CatCatchWgt,
                Depth, meanLat, meanLong, Flag_SweptAreaDSKM2, date,
                Density_kg,Density_kg_log10,HaulDur,SweptAreaDSKM2,CPUE,
                ScientificName_WoRMS,LngtClass,Haul_ID) %>% 
  mutate(Species_group = factor(case_when(LngtClass > 34 & ScientificName_WoRMS == 'Gadus morhua' ~ "Adult Gadus morhua",
                                          LngtClass < 35 & ScientificName_WoRMS == 'Gadus morhua' ~ "Juvenile Gadus morhua",
                                          TRUE ~ ScientificName_WoRMS))) 

# Create trawl level data to model and map trawl level indices
dat_trawl <- dat %>% 
  group_by(Haul_ID, Species_group, date, Year, HaulDur,
           Depth, meanLat, meanLong, Quarter, Gear, .drop = FALSE) %>% 
  summarise(Density_kg = sum(Density_kg),
            SweptAreaDSKM2 = sum(SweptAreaDSKM2),
            CPUE = sum(CPUE),
            HaulDur = sum(HaulDur)) %>% 
  ungroup() %>% 
  fill(Haul_ID, Species_group, date, Year, HaulDur,
       Depth, meanLat, meanLong, Quarter, Gear) %>% 
  arrange(date, meanLat, meanLong) %>% 
  filter(Species_group %in% taxa)  %>% 
  filter(Year !=2008 | Species_group != "Juvenile Gadus morhua") %>% 
  filter(Year !=2008 | Species_group != "Adult Gadus morhua") %>% 
  filter(Year !=2008 | Species_group != "Platichthys flesus") %>% 
  mutate(ScientificName_WoRMS = factor(Species_group, 
                                       levels = c("Limanda limanda",
                                                  "Platichthys flesus",
                                                  "Pleuronectes platessa",
                                                  "Juvenile Gadus morhua",
                                                  "Adult Gadus morhua"),
                                       labels = c("Dab",
                                                  "Flounder",
                                                  "Plaice",
                                                  "Juvenile Cod",
                                                  "Adult Cod"))) 


####################################################################
# Plotting BITS features
####################################################################

TS <- dat_trawl %>% 
  ggplot(aes(Year, log(Density_kg+1), colour = factor(Quarter))) +
  stat_summary(fun.data = "mean_cl_boot",
               alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  guides(color=guide_legend(override.aes = list(size=1.5)))+
  labs(x = "Year", y = expression(log(kg~km^-2)),
       colour = "Quarter")+
  facet_grid(~ScientificName_WoRMS, scales = "free_y")+
  theme_bw(16)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8, size=10),
        strip.text = element_text(size = 10))
ggsave("../Figures/Spring2024Revision/Timeseries_Biomass.png", TS, 
       width = 17, height = 15, units = "cm", dpi = 600)  


indiv_v_weight <- dat_trawl %>% 
  filter(Gear != "NA") %>% 
  ggplot(aes(x = log(CPUE+1), y = log((Density_kg)+1))) +
  geom_point(aes(color = Gear), alpha=0.7, size=2)+
  #geom_smooth(method = "loess")+
  facet_grid(~ScientificName_WoRMS)+
  labs(x= "log(Individuals per hour)", y = expression(log(kg~km^2)))+
  ggpubr::stat_cor(aes(label = ..r.label..),
                   method = "pearson", size=4, label.y = 8.5,
                    r.accuracy = 0.01)+
  scale_color_brewer(palette = "Set1")+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme_bw(16)+
  theme(axis.text.x = element_text(angle = 35, hjust=0.8, size=10),
        strip.text = element_text(size = 10))+
  scale_y_continuous(limits = c(0, 9))
ggsave("../Figures/Spring2024Revision/Weight_v_catch.png", indiv_v_weight, 
       width = 17, height = 15, units = "cm", dpi = 600)  

Quarter <- c("Q1", "Q4")
names(Quarter) <- c(1, 4)
Swept_v_HaulDur <- dat_trawl %>% 
  filter(SweptAreaDSKM2 < 50) %>% 
  ggplot(aes(x = HaulDur, y = SweptAreaDSKM2)) +
  geom_point(aes(color = Gear), alpha=0.7, size=2)+
  geom_smooth(method = "lm")+
  facet_wrap(~Quarter,labeller = labeller(Quarter=Quarter))+
  labs(x= "Haul Duration (minutes)", y = expression(Swept~Area~(km^2)))+
  ggpubr::stat_cor(aes(label = ..r.label..),
                   method = "pearson", size=4, label.x = 2,
                   r.accuracy = 0.01)+
  scale_color_brewer(palette = "Set1")+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme_bw(16)
ggsave("../Figures/Spring2024Revision/Swept_v_Dur.png", Swept_v_HaulDur, 
       width = 17, height = 15, units = "cm", dpi = 600)  


depthswept <- dat_trawl %>% 
  filter(SweptAreaDSKM2 < 50) %>% 
  filter(Depth < 126) %>% 
  ggplot(aes(x = Depth, y = SweptAreaDSKM2)) +
  geom_point(aes(color = Gear), alpha=0.7, size=2)+
  facet_wrap(~Quarter,labeller = labeller(Quarter=Quarter))+
  labs(x= "Depth (m)", y = expression(Swept~Area~(km^2)))+
  scale_color_brewer(palette = "Set1")+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme_bw(16)
ggsave("../Figures/Spring2024Revision/Swept_v_Depth.png", depthswept, 
       width = 17, height = 15, units = "cm", dpi = 600)  

