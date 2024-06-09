library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(ggpubr)

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

Coastline <- read_sf("../../../Feb2023_Transfer/Oceanographic/Data/GSHHS_shp/f/", layer = "GSHHS_f_L1") %>% 
  st_transform(4326) %>% 
  st_crop(st_bbox(icesarea))

# Extra for cropping tiled sf smooths
icesarea_sf <- read_sf("../../Oceanographic/Data/ICES_areas/", layer = "ICES_Areas_20160601_cut_dense_3857") %>%
  filter(SubDivisio == "21" |SubDivisio == "22" | SubDivisio == "23" |SubDivisio == "24"|SubDivisio == "25"|
           SubDivisio == "26"|SubDivisio == "27"|SubDivisio == "28"|SubDivisio == "29") %>% 
  st_transform(st_crs(4326)) %>% 
  st_crop(xmin=xlim[1], ymin=ylim[1], xmax=xlim[2], ymax=ylim[2]) 

################################################################
# Approach 1
################################################################

GAM <- mgcv::bam(Density_log ~ 
                   te(Latitude, Longitude, Year, by=ScientificName_WoRMS) + 
                   #s(Year, Quarter, by = ScientificName_WoRMS,  bs="fs",k=20)+ 
                   te(Year, Depth, by = ScientificName_WoRMS, m=2, bs = "fs")+
                   s(ScientificName_WoRMS, bs = "re"),
                 data = dat, 
                 family = tw(), 
                 method = 'fREML', 
                 select = T,
                 discrete=T,
                 rho = 0.1713523)
#(valRho <- acf(resid(GAM), plot=FALSE)$acf[2])

par(mfrow = c(2,2))
gam.check(GAM)
round(k.check(GAM),3)

par(mfrow = c(2,2))
appraise(GAM)&theme_bw()

summary(GAM)

par(mfrow = c(3,4))
plot(GAM, shade = T)

par(mfrow=c(1,2))
acf(resid(GAM), main = "ACF")
pacf(resid(GAM), main = "pACF")
dev.off()

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
                     #s(Year, Quarter, by = ScientificName_WoRMS,  bs="fs",k=20)+ 
                     te(Year, Depth, by = ScientificName_WoRMS, m=2, bs = "fs")+
                     s(ScientificName_WoRMS, bs = "re"),
                   data = train_data, 
                   family = tw(), 
                   method = 'fREML', 
                   select = T,
                   discrete=T,
                   rho = 0.1713523)
  
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

aic_df <- do.call(rbind, AIC)
mae_df <- do.call(rbind, mae)
disp_df <- do.call(rbind, disp)
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
  filter(FoldID !=10) %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(disp_mean = mean(disp),
            SD_dis = sd(disp))
df_result %>% 
  filter(FoldID !=10) %>% 
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
ggsave("../Figures/Spring2024Revision/Correlation_10fold_Species.png", width = 17,
       height = 15, units = "cm", dpi = 600)


pred_disp <- df_result %>% 
  #full_join(dat) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= disp, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  geom_hline(yintercept =  1, linetype=2, linewidth =1)+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Precision (Dispersion)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/Dispersion_10fold_Species.png", width = 17,
       height = 15, units = "cm", dpi = 600)

pred_mae <- df_result %>% 
  #full_join(dat) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= mae, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  geom_hline(yintercept =  0, linetype=2, linewidth =1)+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Accuracy (MAE)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/MAE_10fold_Species.png", width = 17,
       height = 15, units = "cm", dpi = 600)

############################################################################
# Model- II  geography + abiotic with species-specific smoothers
############################################################################

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
                  Depth,Density_log, FoldID,
                  bottomT, BottomOxygen, BottomSalinity)
  
  train_data <- dat_m %>% dplyr::filter(FoldID != fold)
  val_data <- dat_m %>% dplyr::filter(FoldID == fold)
  
  GAM_2 <- mgcv::bam(Density_log ~ 
                       te(Latitude, Longitude, Year, by=ScientificName_WoRMS) + 
                       s(Year, by = ScientificName_WoRMS)+ 
                       s(Depth, ScientificName_WoRMS, bs="fs", m=2,k=45)+
                       s(bottomT, ScientificName_WoRMS, bs="fs", m=2,k=45)+
                       #s(chl, ScientificName_WoRMS, bs="fs",m=2,k=45)+
                       s(BottomOxygen, ScientificName_WoRMS, bs="fs",k=45,m=2)+
                       s(BottomSalinity, ScientificName_WoRMS, bs="fs",m=2)+
                       s(ScientificName_WoRMS, bs = "re"),
                     data = dat, 
                     family = tw(), 
                     method = 'fREML', 
                     discrete = T,
                     select = T,
                     rho = 0.163319)
  
  pred <- predict(GAM_2, newdata = val_data, type="response")
  
  #CO2_modG_pred <- cbind(val_pred,
  #                       predict(GAM, 
  #                               CO2_modG_pred, 
  #                               se.fit=TRUE, 
  #                               type="response"))
  
  # Calculate metrics 
  #AIC
  AIC <- append(aic, list(AIC(GAM_2)))
  
  # MAE (accuracy)
  #print(mean(abs(dat_test_fold$Density_log10 - pred)))
  mae <- append(mae, list(mean(abs(val_data$Density_log - pred))))
  
  preds <- data.frame(val_data$FoldID,
                      val_data$Density_log, 
                      pred, 
                      val_data$ScientificName_WoRMS,
                      val_data$Year,
                      val_data$Quarter,
                      val_data$bottomT, 
                      val_data$BottomOxygen, 
                      val_data$BottomSalinity)
  
  # Discrimination (Pearson correlation - scaled and unscaled)
  cors <- append(cors, list(preds))
  
  # dispersion (precision) based on Waldock et al. 2022
  disp <- append(disp, list(var(pred)/var(val_data$Density_log)))
  
}

#aic_df <- do.call(rbind, AIC)
#mae_df <- do.call(rbind, mae)
#disp_df <- do.call(rbind, disp)
cors_df_gam2 <- do.call(rbind, cors)

colnames(cors_df_gam2) <- c("FoldID","Density_log",
                            "Predicted", "ScientificName_WoRMS",
                            "Year", "Quarter",
                            "bottomT", "BottomOxygen", "BottomSalinity")

folds <- unique(cors_df_gam2$FoldID)
Species <- unique(cors_df_gam2$ScientificName_WoRMS)
df_ls_gam2 <- list()

for(i in folds){
  for(q in Species){
    x <- filter(cors_df_gam2, FoldID == i & ScientificName_WoRMS == q) %>% 
      mutate(disp = var(Predicted)/var(Density_log),
             mae = mean(abs(Density_log-Predicted)))
    
    df_ls_gam2 <- append(df_ls_gam2, list(x))
  }
}

df_result_gam2 <- do.call(rbind, df_ls_gam2)


#############################################################
#
df_result_gam2 %>% 
  group_by(ScientificName_WoRMS, FoldID) %>% 
  mutate(cors = cor(Density_log, Predicted)) %>% 
  ungroup() %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(cor_mean = mean(cors),
            SD_cor = sd(cors))


df_result_gam2 %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(disp_mean = mean(disp),
            SD_dis = sd(disp))

df_result_gam2 %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(mae_mean = mean(mae),
            SD_mae = sd(mae))

#############################################################

pred_cor <- df_result_gam2 %>% 
  dplyr::mutate(FoldID = factor(FoldID)) %>% 
  #dplyr::mutate(Density_kg = exp(Density_log),
  #       Predicted_kg = exp(Predicted)) %>% 
  ggplot(aes(x = Density_log, y= Predicted))+
  geom_point(aes(color = FoldID), size = 2, alpha = 0.5)+
  geom_abline(slope=1, intercept=0)+
  theme_bw(14)+
  theme(strip.text = element_text(size=10))+
  scale_y_continuous(limits = c(0,8.5))+
  facet_wrap(~ScientificName_WoRMS)+
  labs(x = expression(Observed~log(kg~km^-2)), 
       y = expression(Predicted~log(kg~km^-2)))+
  stat_smooth(aes(color = FoldID), method=lm)+
  stat_cor(aes(color = FoldID, label = ..r.label..),
           method = "pearson", 
           size=2.5,  r.accuracy = 0.01, 
           label.x =  0)+
  scale_color_brewer(palette = "Paired")

ggsave("../Figures/Spring2024Revision/Correlation_GAM2_10fold_Species.png",
       width = 17,height = 15, units = "cm", dpi = 600)


pred_disp <- df_result_gam2 %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  #filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= disp, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  geom_hline(yintercept =  1, linetype=2, linewidth =1)+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Precision (Dispersion)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/Dispersion_GAM2_10fold_Species.png", pred_disp,
       width = 17,height = 15, units = "cm", dpi = 600)

pred_mae <- df_result_gam2 %>% 
  #full_join(dat) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  #filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= mae, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  geom_hline(yintercept =  0, linetype=2, linewidth =1)+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Accuracy (MAE)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/MAE_GAM2_10fold_Species.png", pred_mae, 
       width = 17,height = 15, units = "cm", dpi = 600)


############################################################################
# Model- III  geography + abiotic with seasonal independence
############################################################################

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
                  Depth,Density_log, FoldID,
                  bottomT, BottomOxygen, BottomSalinity)
  
  train_data <- dat_m %>% dplyr::filter(FoldID != fold)
  val_data <- dat_m %>% dplyr::filter(FoldID == fold)
  
  GAM_3 <- mgcv::bam(Density_log ~ 
                       te(Latitude, Longitude, Year, by = ScientificName_WoRMS) + 
                       s(Year, ScientificName_WoRMS, by = Quarter,bs="fs")+ 
                       s(Depth, ScientificName_WoRMS, by = Quarter, bs="fs", m=2,k=45)+
                       s(bottomT, ScientificName_WoRMS, by = Quarter,bs="fs", m=2,k=45)+
                       #s(chl, ScientificName_WoRMS, by = Quarter, bs="fs",m=2,k=45)+
                       s(BottomOxygen, ScientificName_WoRMS, by = Quarter, bs="fs", m=2, k=45)+
                       s(BottomSalinity, ScientificName_WoRMS, by = Quarter, m=2, bs="fs")+
                       s(ScientificName_WoRMS, bs = "re"),
                     data = dat, 
                     family = tw(), 
                     method = 'fREML', 
                     select = T,
                     discrete = T,
                     rho=0.1526755
  )
  
  pred <- predict(GAM_3, newdata = val_data, type="response")
  
  #CO2_modG_pred <- cbind(val_pred,
  #                       predict(GAM, 
  #                               CO2_modG_pred, 
  #                               se.fit=TRUE, 
  #                               type="response"))
  
  # Calculate metrics 
  #AIC
  AIC <- append(aic, list(AIC(GAM_3)))
  
  # MAE (accuracy)
  #print(mean(abs(dat_test_fold$Density_log10 - pred)))
  mae <- append(mae, list(mean(abs(val_data$Density_log - pred))))
  
  preds <- data.frame(val_data$FoldID,
                      val_data$Density_log, 
                      pred, 
                      val_data$ScientificName_WoRMS,
                      val_data$Year,
                      val_data$Quarter,
                      val_data$bottomT, 
                      val_data$BottomOxygen, 
                      val_data$BottomSalinity)
  
  # Discrimination (Pearson correlation - scaled and unscaled)
  cors <- append(cors, list(preds))
  
  # dispersion (precision) based on Waldock et al. 2022
  disp <- append(disp, list(var(pred)/var(val_data$Density_log)))
  
}

#aic_df <- do.call(rbind, AIC)
#mae_df <- do.call(rbind, mae)
#disp_df <- do.call(rbind, disp)
cors_df_gam3 <- do.call(rbind, cors)

colnames(cors_df_gam3) <- c("FoldID","Density_log",
                            "Predicted", "ScientificName_WoRMS",
                            "Year", "Quarter",
                            "bottomT", "BottomOxygen", "BottomSalinity")

folds <- unique(cors_df_gam3$FoldID)
Species <- unique(cors_df_gam3$ScientificName_WoRMS)
df_ls_gam3 <- list()

for(i in folds){
  for(q in Species){
    x <- filter(cors_df_gam3, FoldID == i & ScientificName_WoRMS == q) %>% 
      mutate(disp = var(Predicted)/var(Density_log),
             mae = mean(abs(Density_log-Predicted)))
    
    df_ls_gam3 <- append(df_ls_gam3, list(x))
  }
}

df_result_gam3 <- do.call(rbind, df_ls_gam3)


#############################################################
#
df_result_gam3 %>% 
  group_by(ScientificName_WoRMS, FoldID) %>% 
  mutate(cors = cor(Density_log, Predicted)) %>% 
  ungroup() %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(cor_mean = mean(cors),
            SD_cor = sd(cors))

  
df_result_gam3 %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(disp_mean = mean(disp),
            SD_dis = sd(disp))
df_result_gam3 %>% 
  group_by(ScientificName_WoRMS) %>% 
  summarise(mae_mean = mean(mae),
            SD_mae = sd(mae))

#############################################################

pred_cor <- df_result_gam3 %>% 
  dplyr::mutate(FoldID = factor(FoldID)) %>% 
  #dplyr::mutate(Density_kg = exp(Density_log),
  #       Predicted_kg = exp(Predicted)) %>% 
  ggplot(aes(x = Density_log, y= Predicted))+
  geom_point(aes(color = FoldID), size = 2, alpha = 0.5)+
  geom_abline(slope=1, intercept=0)+
  theme_bw(14)+
  theme(strip.text = element_text(size=10))+
  scale_y_continuous(limits = c(0,8.5))+
  facet_wrap(~ScientificName_WoRMS)+
  labs(x = expression(Observed~log(kg~km^-2)), 
       y = expression(Predicted~log(kg~km^-2)))+
  stat_smooth(aes(color = FoldID), method=lm)+
  stat_cor(aes(color = FoldID, label = ..r.label..),
           method = "pearson", 
           size=2.5,  r.accuracy = 0.01, 
           label.x =  0)+
  scale_color_brewer(palette = "Paired")

ggsave("../Figures/Spring2024Revision/Correlation_GAM3_10fold_Species.png",
       width = 17,height = 15, units = "cm", dpi = 600)


pred_disp <- df_result_gam3 %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  #filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= disp, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  #facet_wrap(~ScientificName_WoRMS)+
  geom_hline(yintercept =  1, linetype=2, linewidth =1)+
  labs(y = "Precision (Dispersion)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/Dispersion_GAM3_10fold_Species.png", pred_disp,
       width = 17,height = 15, units = "cm", dpi = 600)

pred_mae <- df_result_gam3 %>% 
  #full_join(dat) %>% 
  mutate(FoldID = factor(FoldID)) %>% 
  #filter(FoldID != 10) %>%
  ggplot(aes(x = (ScientificName_WoRMS), y= mae, fill = ScientificName_WoRMS))+
  geom_boxplot()+
  theme_bw(14)+
  geom_hline(yintercept =  0, linetype=2, linewidth =1)+
  #facet_wrap(~ScientificName_WoRMS)+
  labs(y = "Accuracy (MAE)", x="")+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust=1))
ggsave("../Figures/Spring2024Revision/MAE_GAM3_10fold_Species.png", pred_mae, 
       width = 17,height = 15, units = "cm", dpi = 600)
