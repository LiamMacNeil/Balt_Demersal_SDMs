library(tidyverse)
library(vroom)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#######################################################
# All performance metrics  
#######################################################
BITS_Taxa <- c("Flounder", "CodAdult","AdultCod", "CodJuvenile", "JuvenileCod", "Dab", "Plaice")
BITS_Taxa_x <- c("Flounder", "Adult Cod", "Juvenile Cod", "Dab", "Plaice")

gam_metrics <- (list.files("../Model_Metrics/GAM", pattern = "\\.csv$", full.names = T))
rf_metrics <- (list.files("../Model_Metrics/RF", pattern = "\\.csv$", full.names = T))
brt_metrics <- (list.files("../Model_Metrics/BRT", pattern = "\\.csv$", full.names = T))

gam_metrics <- vroom(gam_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa_x[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa_x[4],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa_x[5]))

gam_metrics$Common_name[gam_metrics$Common_name == "AdultCod"] <- "CodAdult"
gam_metrics$Common_name[gam_metrics$Common_name == "JuvenileCod"] <- "CodJuvenile"

rf_metrics <- vroom(rf_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa_x[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa_x[4],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa_x[5])) 

rf_metrics$Common_name[rf_metrics$Common_name == "JuvenileCod"] <- "CodJuvenile"

# Clean vroom method robust to new data entries
brt_metrics <- vroom(brt_metrics, id = 'filename') %>% 
  mutate(Common_name = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa_x[1],
                                 grepl(BITS_Taxa[2], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[3], filename) ~ BITS_Taxa_x[2],
                                 grepl(BITS_Taxa[4], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[5], filename) ~ BITS_Taxa_x[3],
                                 grepl(BITS_Taxa[6], filename) ~ BITS_Taxa_x[4],
                                 grepl(BITS_Taxa[7], filename) ~ BITS_Taxa_x[5])) 


metrics <-full_join(brt_metrics, rf_metrics)
metrics <-full_join(metrics, gam_metrics)

#metrics[metrics == "GLM"] <- "GAM"
metrics$model <- factor(metrics$model, levels = c("GAM", "RF", "BRT"))
metrics$Common_name <- factor(metrics$Common_name, levels = c("Dab",
                                                              "Flounder", 
                                                              "Plaice",
                                                              "Juvenile Cod",
                                                              "Adult Cod"))
metrics <- metrics %>% 
  mutate(Quarter_factor = case_when(Quarter == 1 ~ "Q1",
                                    Quarter == 4 ~ "Q4"))
# Vector of facet labels
Quarters <- c("Q1", "Q4", "2001-2010", "2011-2020")

discrimination <- ggplot(metrics, aes(Common_name, discrimination, fill=model)) +
  geom_boxplot(width = 0.8)+
  facet_grid(Year~Quarter_factor, scales = "free_x")+
  theme_bw(16)+
  #labs(x = "", y = expression(Discrimination~(Pearson*paste("'s")~italic(r))))+
  labs(x = "", y = "")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill = guide_legend("Model"))
ggsave("../Figures/Discrimination.png", discrimination, width = 16, height = 12, units = "cm", dpi = 600)

metrics %>% 
  filter(dispersion < 2) %>%
  group_by(Quarter) %>% 
  summarise(mean(discrimination), mean(mae), mean(dispersion))

MAE <- ggplot(metrics, aes(Common_name, mae, fill=model)) +
  geom_boxplot(width = 0.8)+
  facet_grid(Year~Quarter_factor, scales = "free_x")+
  theme_bw(16)+
  #labs(x = "", y = "Accuracy (MAE)")+
  labs(x = "", y = "")+
  #scale_fill_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 35, hjust =1),
        strip.text = element_text(size = 14))+
  guides(fill = guide_legend("Model"))
ggsave("../Figures/MAE.png", MAE, width = 16, height = 12, units = "cm", dpi = 600)

Dispersion <- metrics %>% filter(dispersion < 2) %>% 
  ggplot(aes(Common_name, dispersion, fill=model)) +
  geom_boxplot(width = 0.8)+
  facet_grid(Year~Quarter_factor, scales = "free_x")+
  theme_bw(16)+
  #labs(x = "", y = "Precision (Dispersion)")+
  labs(x = "", y = "")+
  scale_fill_brewer(palette = "Set1")+
  #scale_fill_brewer(palette = "Paired")+
  theme(axis.text.x = element_text(angle = 35, hjust =1),
        strip.text = element_text(size = 14))+
  guides(fill = guide_legend("Model"))
ggsave("../Figures/Dispersion.png", Dispersion, width = 16, height = 12, units = "cm", dpi = 600)
