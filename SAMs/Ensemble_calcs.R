library(terra)
library(raster)
library(vroom)
library(tidyverse)
# Ensemble approach

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Weighting by pearson correlation coefficient

# Weighted ensemble = (Cor X pred)/Cor (for site i and model j)

# Weighted raster maps

 
#######################################################
# All performance metrics  
#######################################################
BITS_Taxa <- c("Flounder", "CodAdult","AdultCod", "CodJuvenile", "JuvenileCod", "Dab", "Plaice")

gam_metrics <- (list.files("../Model_Metrics/GAM", pattern = "\\.csv$", full.names = T))
rf_metrics <- (list.files("../Model_Metrics/RF", pattern = "\\.csv$", full.names = T))
brt_metrics <- (list.files("../Model_Metrics/BRT", pattern = "\\.csv$", full.names = T))

gam_metrics <- vroom(gam_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa[3],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa[4],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa[5],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa[6],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa[7]))
          
gam_metrics$Common_name[gam_metrics$Common_name == "AdultCod"] <- "CodAdult"
gam_metrics$Common_name[gam_metrics$Common_name == "JuvenileCod"] <- "CodJuvenile"

rf_metrics <- vroom(rf_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa[3],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa[4],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa[5],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa[6],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa[7])) 

# Clean vroom method robust to new data entries
brt_metrics <- vroom(brt_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa[3],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa[4],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa[5],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa[6],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa[7])) 
  


metrics <-full_join(brt_metrics, rf_metrics)
metrics <-full_join(metrics, gam_metrics)

#######################################################
# All Prediction Maps  
#######################################################

layers <- rast(list.files("../Predictions/", 
                           pattern = "\\.tif$", full.names = T))

# Fixing naming scheme
names <- layers %>% 
  sources() %>%
  as.data.frame() %>%
  rename(Column = colnames(layers %>% 
                            sources() %>%
                            as.data.frame())) %>% 
  separate_wider_delim(Column, delim = "/", names = c(".", "Dir",
                                                                                        "Drive", 
                                                                                        "Folder", 
                                                                                        "Models", 
                                                                                        "Files", 
                                                                                        "Layer"), 
                                                         too_few  = "align_end", too_many = "merge")

names(layers) <- c(names$Layer)
#######################################################
# 2001 - 2010
#######################################################

###
# Q1
###
layers[[varnames(layers) == "CodAdult_Q1_2001_GAM"]]

AdultCod_stacked_Q1_2001 <- c(layers[[varnames(layers) == "CodAdult_Q1_2001_GAM"]] ,
                               layers[[varnames(layers) == "CodAdult_Q1_2001_RF"]],
                                   layers[[varnames(layers) == "CodAdult_Q1_2001_BRT"]]
                               )


AdultCod_weighted_Q1_2001 <- weighted.mean(AdultCod_stacked_Q1_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                            dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                            dplyr::filter(model == "RF") %>% 
                                                                            dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                            summarise(median = median(discrimination)) ,
                                                                          metrics %>%  #BRT
                                                                            dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                            dplyr::filter(model == "BRT") %>% 
                                                                            dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                            summarise(median = median(discrimination)) ,
                                                                          metrics %>%  #GAM
                                                                            dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                            dplyr::filter(model == "GAM") %>% 
                                                                            dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                            summarise(median = median(discrimination))))))

JuvenileCod_stacked_Q1_2001 <- c( layers[[varnames(layers) == "CodJuvenile_Q1_2001_RF"]],
                               layers[[varnames(layers) == "CodJuvenile_Q1_2001_BRT"]],
                               layers[[varnames(layers) == "CodJuvenile_Q1_2001_GAM"]])


JuvenileCod_weighted_Q1_2001 <- weighted.mean(JuvenileCod_stacked_Q1_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                                  dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                  dplyr::filter(model == "RF") %>% 
                                                                                                                  dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                                  summarise(median = median(discrimination)) ,
                                                                                                                metrics %>%  #BRT
                                                                                                                  dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                  dplyr::filter(model == "BRT") %>% 
                                                                                                                  dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                                  summarise(median = median(discrimination)),
                                                                                                               metrics %>%  #GAM
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                                 summarise(median = median(discrimination))))))


Plaice_stacked_Q1_2001 <- c( layers[[varnames(layers) == "Plaice_Q1_2001_RF"]],
                               layers[[varnames(layers) == "Plaice_Q1_2001_BRT"]],
                             layers[[varnames(layers) == "Plaice_Q1_2001_GAM"]])


Plaice_weighted_Q1_2001 <- weighted.mean(Plaice_stacked_Q1_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "GAM") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)) ))))


Flounder_stacked_Q1_2001 <- c( layers[[varnames(layers) == "Flounder_Q1_2001_RF"]],
                             layers[[varnames(layers) == "Flounder_Q1_2001_BRT"]],
                             layers[[varnames(layers) == "Flounder_Q1_2001_GAM"]])


Flounder_weighted_Q1_2001 <- weighted.mean(Flounder_stacked_Q1_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                       dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "RF") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)) ,
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "BRT") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)),
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "GAM") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)) ))))


Dab_stacked_Q1_2001 <- c( layers[[varnames(layers) == "Dab_Q1_2001_RF"]],
                               layers[[varnames(layers) == "Dab_Q1_2001_BRT"]],
                          layers[[varnames(layers) == "Dab_Q1_2001_GAM"]])


Dab_weighted_Q1_2001 <- weighted.mean(Dab_stacked_Q1_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 1) %>%
                                                                                                 summarise(median = median(discrimination)) ))))


###
# Q4
###

AdultCod_stacked_Q4_2001 <- c( layers[[varnames(layers) == "CodAdult_Q4_2001_RF"]],
                               layers[[varnames(layers) == "CodAdult_Q4_2001_BRT"]],
                               layers[[varnames(layers) == "CodAdult_Q4_2001_GAM"]])


AdultCod_weighted_Q4_2001 <- weighted.mean(AdultCod_stacked_Q4_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ))))

JuvenileCod_stacked_Q4_2001 <- c( layers[[varnames(layers) == "CodJuvenile_Q4_2001_RF"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q4_2001_BRT"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q4_2001_GAM"]])


JuvenileCod_weighted_Q4_2001 <- weighted.mean(JuvenileCod_stacked_Q4_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)) ,
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)),
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)) ))))


Plaice_stacked_Q4_2001 <- c( layers[[varnames(layers) == "Plaice_Q4_2001_RF"]],
                             layers[[varnames(layers) == "Plaice_Q4_2001_BRT"]],
                             layers[[varnames(layers) == "Plaice_Q4_2001_GAM"]])


Plaice_weighted_Q4_2001 <- weighted.mean(Plaice_stacked_Q4_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "RF") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)) ,
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "BRT") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)),
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "GAM") %>% 
                                                                                                       dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)) ))))


Flounder_stacked_Q4_2001 <- c( layers[[varnames(layers) == "Flounder_Q4_2001_RF"]],
                               layers[[varnames(layers) == "Flounder_Q4_2001_BRT"]],
                               layers[[varnames(layers) == "Flounder_Q4_2001_GAM"]])


Flounder_weighted_Q4_2001 <- weighted.mean(Flounder_stacked_Q4_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ))))


Dab_stacked_Q4_2001 <- c( layers[[varnames(layers) == "Dab_Q4_2001_RF"]],
                          layers[[varnames(layers) == "Dab_Q4_2001_BRT"]],
                          layers[[varnames(layers) == "Dab_Q4_2001_GAM"]])


Dab_weighted_Q4_2001 <- weighted.mean(Dab_stacked_Q4_2001, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)) ,
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)),
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                 dplyr::filter(Year == "2001-2010" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)) ))))



#######################################################
# 2011 - 2020
#######################################################

###
# Q1
###
AdultCod_stacked_Q1_2020 <- c( layers[[varnames(layers) == "CodAdult_Q1_2020_RF"]],
                               layers[[varnames(layers) == "CodAdult_Q1_2020_BRT"]],
                               layers[[varnames(layers) == "CodAdult_Q1_2020_GAM"]])


AdultCod_weighted_Q1_2020 <- weighted.mean(AdultCod_stacked_Q1_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ))))

JuvenileCod_stacked_Q1_2020 <- c( layers[[varnames(layers) == "CodJuvenile_Q1_2020_RF"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q1_2020_BRT"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q1_2020_GAM"]])


JuvenileCod_weighted_Q1_2020 <- weighted.mean(JuvenileCod_stacked_Q1_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                                 summarise(median = median(discrimination)) ,
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                                 summarise(median = median(discrimination)),
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                                 summarise(median = median(discrimination)) ))))


Plaice_stacked_Q1_2020 <- c( layers[[varnames(layers) == "Plaice_Q1_2020_RF"]],
                             layers[[varnames(layers) == "Plaice_Q1_2020_BRT"]],
                             layers[[varnames(layers) == "Plaice_Q1_2020_GAM"]])


Plaice_weighted_Q1_2020 <- weighted.mean(Plaice_stacked_Q1_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "RF") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)) ,
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "BRT") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)),
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "GAM") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                       summarise(median = median(discrimination)) ))))


Flounder_stacked_Q1_2020 <- c( layers[[varnames(layers) == "Flounder_Q1_2020_RF"]],
                               layers[[varnames(layers) == "Flounder_Q1_2020_BRT"]],
                               layers[[varnames(layers) == "Flounder_Q1_2020_GAM"]])


Flounder_weighted_Q1_2020 <- weighted.mean(Flounder_stacked_Q1_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                           summarise(median = median(discrimination)) ))))


Dab_stacked_Q1_2020 <- c( layers[[varnames(layers) == "Dab_Q1_2020_RF"]],
                          layers[[varnames(layers) == "Dab_Q1_2020_BRT"]],
                          layers[[varnames(layers) == "Dab_Q1_2020_GAM"]])


Dab_weighted_Q1_2020 <- weighted.mean(Dab_stacked_Q1_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                 summarise(median = median(discrimination)) ,
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                 summarise(median = median(discrimination)) ,
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 1) %>%
                                                                                                 summarise(median = median(discrimination)) ))))


###
# Q4
###

AdultCod_stacked_Q4_2020 <- c( layers[[varnames(layers) == "CodAdult_Q4_2020_RF"]],
                               layers[[varnames(layers) == "CodAdult_Q4_2020_BRT"]],
                               layers[[varnames(layers) == "CodAdult_Q4_2020_GAM"]])


AdultCod_weighted_Q4_2020 <- weighted.mean(AdultCod_stacked_Q4_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("CodAdult",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ))))

JuvenileCod_stacked_Q4_2020 <- c( layers[[varnames(layers) == "CodJuvenile_Q4_2020_RF"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q4_2020_BRT"]],
                                  layers[[varnames(layers) == "CodJuvenile_Q4_2020_GAM"]])


JuvenileCod_weighted_Q4_2020 <- weighted.mean(JuvenileCod_stacked_Q4_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)) ,
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)),
                                                                                                               metrics %>%  #BRT
                                                                                                                 dplyr::filter(grepl("CodJuvenile",Common_name)) %>% 
                                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                                 summarise(median = median(discrimination)) ))))


Plaice_stacked_Q4_2020 <- c( layers[[varnames(layers) == "Plaice_Q4_2020_RF"]],
                             layers[[varnames(layers) == "Plaice_Q4_2020_BRT"]],
                             layers[[varnames(layers) == "Plaice_Q4_2020_GAM"]])


Plaice_weighted_Q4_2020 <- weighted.mean(Plaice_stacked_Q4_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "RF") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)) ,
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "BRT") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)),
                                                                                                     metrics %>%  #BRT
                                                                                                       dplyr::filter(grepl("Plaice",Common_name)) %>% 
                                                                                                       dplyr::filter(model == "GAM") %>% 
                                                                                                       dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                       summarise(median = median(discrimination)) ))))


Flounder_stacked_Q4_2020 <- c( layers[[varnames(layers) == "Flounder_Q4_2020_RF"]],
                               layers[[varnames(layers) == "Flounder_Q4_2020_BRT"]],
                               layers[[varnames(layers) == "Flounder_Q4_2020_GAM"]])


Flounder_weighted_Q4_2020 <- weighted.mean(Flounder_stacked_Q4_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "RF") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ,
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "BRT") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)),
                                                                                                         metrics %>%  #BRT
                                                                                                           dplyr::filter(grepl("Flounder",Common_name)) %>% 
                                                                                                           dplyr::filter(model == "GAM") %>% 
                                                                                                           dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                           summarise(median = median(discrimination)) ))))


Dab_stacked_Q4_2020 <- c( layers[[varnames(layers) == "Dab_Q4_2020_RF"]],
                          layers[[varnames(layers) == "Dab_Q4_2020_BRT"]],
                          layers[[varnames(layers) == "Dab_Q4_2020_GAM"]])


Dab_weighted_Q4_2020 <- weighted.mean(Dab_stacked_Q4_2020, w = as.numeric(as.vector(data.frame(metrics %>%  #RF
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "RF") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)) ,
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "BRT") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)),
                                                                                               metrics %>%  #BRT
                                                                                                 dplyr::filter(grepl("Dab",Common_name)) %>% 
                                                                                                 dplyr::filter(model == "GAM") %>% 
                                                                                                 dplyr::filter(Year == "2011-2020" | Quarter == 4) %>%
                                                                                                 summarise(median = median(discrimination)) ))))




############ ############ ############ 
# Ensemble
############ ############ ############

ensembles <- c(AdultCod_weighted_Q1_2001,
               AdultCod_weighted_Q4_2001,
               AdultCod_weighted_Q1_2020,
               AdultCod_weighted_Q4_2020,
               
               JuvenileCod_weighted_Q1_2001,
               JuvenileCod_weighted_Q4_2001,
               JuvenileCod_weighted_Q1_2020,
               JuvenileCod_weighted_Q4_2020,
               
               Plaice_weighted_Q1_2001,
               Plaice_weighted_Q4_2001,
               Plaice_weighted_Q1_2020,
               Plaice_weighted_Q4_2020,
               
               Flounder_weighted_Q1_2001,
               Flounder_weighted_Q4_2001,
               Flounder_weighted_Q1_2020,
               Flounder_weighted_Q4_2020,
               
               Dab_weighted_Q1_2001,
               Dab_weighted_Q4_2001,
               Dab_weighted_Q1_2020,
               Dab_weighted_Q4_2020
                   )

names(ensembles) <- c("AdultCod_Q1_2001",
                      "AdultCod_Q4_2001",
                      "AdultCod_Q1_2020",
                      "AdultCod_Q4_2020",
                      
                      "JuvenileCod_Q1_2001",
                      "JuvenileCod_Q4_2001",
                      "JuvenileCod_Q1_2020",
                      "JuvenileCod_Q4_2020",
                      
                      "Plaice_Q1_2001",
                      "Plaice_Q4_2001",
                      "Plaice_Q1_2020",
                      "Plaice_Q4_2020",
                      
                      "Flounder_Q1_2001",
                      "Flounder_Q4_2001",
                      "Flounder_Q1_2020",
                      "Flounder_Q4_2020",
                      
                      
                      "Dab_Q1_2001",
                      "Dab_Q4_2001",
                      "Dab_Q1_2020",
                      "Dab_Q4_2020")

writeRaster(ensembles, "../Predictions/Ensemble_BRT_RF_GAM.tif", overwrite = T)
