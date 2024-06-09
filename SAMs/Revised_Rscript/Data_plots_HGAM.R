library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(raster)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#data
dat <- read_csv("../Data/Taxa_env_GAMs_v2.csv") %>% 
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
  filter(Year !=2008 | Density_log < 7) %>% 
  #filter(Year != 2008 & 
  #         ScientificName_WoRMS != "Juvenile Gadus morhua" | 
  #         ScientificName_WoRMS != "Adult Gadus morhua") %>% 
  mutate(Year_fac = factor(Year, levels =c(2001:2020), ordered = T))


ggplot(dat, aes(Density_log, Depth, color = ScientificName_WoRMS))+
  geom_point()+
  scale_y_reverse()+
  theme_bw(16)+
  facet_grid(Year~Quarter)
