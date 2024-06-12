library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(raster)
library(viridis)
library(itsadug)
library(oce)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#data
dat <- read_csv("../Data/Taxa_env_GAMs_v2_cutcovs.csv") %>% 
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
  mutate(Year_fac = factor(Year, levels =c(2001:2020), ordered = T)) %>% 
  as.data.frame()

xlim <- c(9.420922, 22.51811)
ylim <- c(53.91644, 59.23877)

#To be deal some stubbornly invalid geometries in the shapefiles 
sf_use_s2(F)

# High res polygons
icesarea <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) %>% 
  as("Spatial")  

icesarea_sf <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

Coastline <- read_sf("../../Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
  st_transform(4326) %>% 
  st_crop(st_bbox(icesarea))

#icesarea_sf <- read_sf("../../../Oceanographic/Data/ICES_shapefiles/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
#  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
#           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
#  st_transform(st_crs(4326)) %>% 
#  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

# Extra for cropping tiled sf smooths
###########################
BITS_buffer <- dat %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = st_crs(4326)) %>% 
  st_buffer(dist = 0.5) %>% 
  st_geometry() %>% 
  st_union()

#####################################################################
# Model -I Geography-only
#####################################################################
k <- 10
FoldID <- rep(1:k, ceiling(nrow(dat)/k))
length(FoldID) <- nrow(dat)
FoldID <- sample(FoldID, replace = FALSE)

results <- list()

mae <- list()
cors <- list()
disp <- list()
aic <- list()

for(fold in 1:k){
  
  dat_m <- dat %>% 
    mutate(FoldID = FoldID) %>% 
    dplyr::select(Latitude, Longitude, Quarter,
                  Year, ScientificName_WoRMS, 
                  Depth,Density_log, FoldID)
  
  train_data <- dat_m %>% dplyr::filter(FoldID != fold)
  val_data <- dat_m %>% dplyr::filter(FoldID == fold)
  
  GAM <- mgcv::bam(Density_log ~ 
                     te(Latitude, Longitude, Year, by=ScientificName_WoRMS) + 
                     te(Year, Depth, by = ScientificName_WoRMS, m=2, bs = "fs")+
                     s(ScientificName_WoRMS, bs = "re"),
                   data = train_data, 
                   family = tw(), 
                   method = 'fREML', 
                   select = T,
                   discrete=T)
  
  pred <- predict(GAM, newdata = val_data, type="response")
    
  #CO2_modG_pred <- cbind(val_pred,
  #                       predict(GAM, 
  #                               CO2_modG_pred, 
  #                               se.fit=TRUE, 
  #                               type="response"))
  
  # Calculate metrics 
  #AIC
  AIC <- append(aic, list(AIC(GAM)))
  
  # MAE (accuracy)
  #print(mean(abs(dat_test_fold$Density_log10 - pred)))
  mae <- append(mae, list(mean(abs(val_data$Density_log - pred))))
  
  preds <- data.frame(val_data$FoldID,
                      val_data$Density_log, 
                      pred, 
                      val_data$ScientificName_WoRMS,
                      val_data$Year,
                      val_data$Quarter)
  
  # Discrimination (Pearson correlation - scaled and unscaled)
  cors <- append(cors, list(preds))
  
  # dispersion (precision) based on Waldock et al. 2022
  disp <- append(disp, list(var(pred)/var(val_data$Density_log)))
  
}

#aic_df <- do.call(rbind, AIC)
#mae_df <- do.call(rbind, mae)
#disp_df <- do.call(rbind, disp)
cors_df <- do.call(rbind, cors)

colnames(cors_df) <- c("FoldID","Density_log", 
                       "Predicted", "ScientificName_WoRMS",
                       "Year",
                       "Quarter")

folds <- unique(cors_df$FoldID)
Species <- unique(cors_df$ScientificName_WoRMS)
df_ls <- list()

for(i in folds){
  for(q in Species){
    x <- filter(cors_df, FoldID == i & ScientificName_WoRMS == q) %>% 
      mutate(disp = var(Predicted)/var(Density_log),
             mae = mean(abs(Density_log-Predicted)))

    df_ls <- append(df_ls, list(x))
  }
}

df_result <- do.call(rbind, df_ls)

#############################################################
#
df_result %>% 
  filter(FoldID !=10) %>% 
  group_by(ScientificName_WoRMS, FoldID) %>% 
  mutate(cors = cor(Density_log, Predicted)) %>% 
  ungroup() %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(cor_mean = mean(cors),
            SD_cor = sd(cors))

df_result %>% 
  filter(FoldID !=8) %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(disp_mean = mean(disp),
            SD_dis = sd(disp))
df_result %>% 
  filter(FoldID !=8) %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(mae_mean = mean(mae),
            SD_mae = sd(mae))

#############################################################

pred_cor <- df_result %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>% 
  #filter(Density_log <5) %>% 
  mutate(Density_kg = exp(Density_log),
         Predicted_kg = exp(Predicted)) %>% 
  ggplot(aes(x = (Density_log), y= Predicted))+
  geom_point(aes(color = FoldID), size = 2, alpha = 0.5)+
  geom_abline(slope=1, intercept=0)+
  theme_bw(14)+
  theme(strip.text = element_text(size=10))+
  facet_wrap(~ScientificName_WoRMS)+
  labs(x = expression(Observed~log(kg~km^-2)), y = expression(Predicted~log(kg~km^-2)))+
  stat_smooth(aes(color = FoldID), method=lm)+
  stat_cor(aes(color = FoldID, label = ..r.label..),
           method = "pearson", 
           size=2.5,  r.accuracy = 0.01, 
           label.x =  0)+
  scale_color_brewer(palette = "Paired")
ggsave("../Figures/Spring2024Revision/Correlation_10fold_Species.png", pred_cor,
       width = 17,height = 15, units = "cm", dpi = 600)

###########################################################################
pred_cor_box <- df_result %>% 
  mutate(Group = case_when(ScientificName_WoRMS == "Pleuronectes platessa" ~ "Plaice",
                           ScientificName_WoRMS == "Limanda limanda" ~ "Dab",
                           ScientificName_WoRMS == "Platichthys flesus" ~ "Flounder",
                           ScientificName_WoRMS == "Juvenile Gadus morhua" ~ "Juvenile Cod",
                           ScientificName_WoRMS == "Adult Gadus morhua" ~ "Adult Cod",
                           TRUE ~ ScientificName_WoRMS)) %>% 
  mutate(Group = factor(Group, levels= c("Dab", "Flounder", "Plaice",
                                         "Juvenile Cod", "Adult Cod"))) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>%
  group_by(ScientificName_WoRMS, FoldID) %>% 
  mutate(cor = cor(Density_log, Predicted)) %>% 
  ggplot(aes(x = (Group), y= cor, fill = Group))+
  geom_boxplot()+
  theme_bw(12)+
  geom_hline(yintercept =  1, linetype=2, linewidth =1)+
  scale_y_continuous(limits = c(0,1))+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Discrimination (Correlation)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust=1))
ggsave("../Figures/Spring2024Revision/CorrelationBox_10fold_Species.png", pred_cor_box,
       width = 17,height = 15, units = "cm", dpi = 600)


pred_disp <- df_result %>% 
  mutate(Group = case_when(ScientificName_WoRMS == "Pleuronectes platessa" ~ "Plaice",
                           ScientificName_WoRMS == "Limanda limanda" ~ "Dab",
                           ScientificName_WoRMS == "Platichthys flesus" ~ "Flounder",
                           ScientificName_WoRMS == "Juvenile Gadus morhua" ~ "Juvenile Cod",
                           ScientificName_WoRMS == "Adult Gadus morhua" ~ "Adult Cod",
                           TRUE ~ ScientificName_WoRMS)) %>% 
  mutate(Group = factor(Group, levels= c("Dab", "Flounder", "Plaice",
                                         "Juvenile Cod", "Adult Cod"))) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>%
  ggplot(aes(x = (Group), y= disp, 
             fill = Group))+
  geom_boxplot()+
  theme_bw(12)+
  geom_hline(yintercept =  1, linetype=2, linewidth =1)+
  scale_y_continuous(limits = c(0,1))+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Precision (Dispersion)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust=1))
ggsave("../Figures/Spring2024Revision/Dispersion_10fold_Species.png", pred_disp,
       width = 17,height = 15, units = "cm", dpi = 600)

pred_mae <- df_result %>% 
  mutate(Group = case_when(ScientificName_WoRMS == "Pleuronectes platessa" ~ "Plaice",
                           ScientificName_WoRMS == "Limanda limanda" ~ "Dab",
                           ScientificName_WoRMS == "Platichthys flesus" ~ "Flounder",
                           ScientificName_WoRMS == "Juvenile Gadus morhua" ~ "Juvenile Cod",
                           ScientificName_WoRMS == "Adult Gadus morhua" ~ "Adult Cod",
                           TRUE ~ ScientificName_WoRMS)) %>% 
  mutate(Group = factor(Group, levels= c("Dab", "Flounder", "Plaice",
                                         "Juvenile Cod", "Adult Cod"))) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>%
  ggplot(aes(x = (Group), y= mae, 
             fill = Group))+
  geom_boxplot()+
  theme_bw(12)+
  geom_hline(yintercept =  0, linetype=2, linewidth =1)+
  scale_y_continuous(limits = c(0,1))+
  labs(y = "Accuracy (MAE)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust=1))
ggsave("../Figures/Spring2024Revision/MAE_10fold_Species.png", pred_mae,
       width = 17,height = 15, units = "cm", dpi = 600)


patch <- pred_cor_box + pred_disp + pred_mae + plot_layout(ncol = 3, nrow = 1)
ggsave("../Figures/Spring2024Revision/Performance_Species.png", patch,
       width = 17,height = 15, units = "cm", dpi = 600)

