library(tidyverse)
library(lubridate)
library(mapdata)
library(gisland)
library(mgcv)
library(ggeffects)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

totalwgt_v <- read_csv("TotalWgt_Valid_and_Pelagic.csv")

# Missing depth values?
totalwgt_v_1 <- totalwgt_v %>% filter(Depth != "NA") %>% 
  filter(CatchWgt != "NA")

# Extreme haul duration values?
totalwgt_v_1 <- totalwgt_v_1 %>% 
  filter(HaulDur >! 90) %>% 
  filter(HaulDur >= 15)

# Calculating distance manually? 
totalwgt_v_1_dist <- totalwgt_v_1 %>% 
  filter(Distance !="NA") %>% 
  filter(Taxa == "Cod") %>% 
  filter(Year > 2001)

totalwgt_v_1_dist <- totalwgt_v_1_dist %>% 
  dplyr::mutate(perarea = ifelse(DataType=="R" | DataType == "S",
                                 CatchWgt / Distance,
                                 CatchWgt))

model_lm <- gam(CatchWgt~Year+ Gear + Depth, data = totalwgt_v_1_dist)

# smoothing parameter with cubmic regression splines
model_gam1 <- gam(CatchWgt ~ s(Year) + s(Depth) + Gear, data = totalwgt_v_1_dist)

AIC(model_lm)
AIC(model_gam1)

summary(model_lm)
summary(model_gam1)

anova(model_lm, model_gam1, test = "Chisq")

plot(ggeffects::ggpredict(model_gam1, facets = T))
vis.gam(model_gam1, type = 'response', plot.type = 'contour')

model_gam2 <- gam(CatchWgt ~ te(Year, Depth, k = 4) , data = totalwgt_v_1_dist)

summary(model_gam2)
plot(ggeffects::ggpredict(model_gam2, facets = T))
vis.gam(
     model_gam2,
     type      = 'response',
     plot.type = 'persp',
     phi       = 30,
     theta     = 30,
     n.grid    = 500,
     border    = NA
   )

## GAMM

model_gamm1 <- bam(perarea ~ te(Year, Depth, Distance, bs = "cs")+ Month_name, random = Gear , data = totalwgt_v_1_dist)

