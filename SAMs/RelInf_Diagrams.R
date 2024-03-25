library(vroom)
library(tidyverse)
library(stringr)
library(oce)
#remotes::install_github("juliasilge/tidytext")
library(tidytext)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

relinfs <- (list.files("../Data/BITS_Var_RelInf", pattern = "\\.csv$", full.names = T))

BITS_Taxa <- c("Flounder", "Cod", "Dab", "Plaice")

# Clean vroom method robust to new data entries
relinfs <- vroom(relinfs,id = 'filename', col_select = c(filename, var, rel.inf)) %>% 
  mutate(Decade = ifelse(str_detect(filename, "2001"), "2001-2010", "2011-2020")) %>% 
  mutate(Quarter = ifelse(str_detect(filename, "Q1"), "Q1", "Q4")) %>% 
  mutate(Taxa = case_when(grepl(BITS_Taxa[1], filename) ~ BITS_Taxa[1],
                          grepl(BITS_Taxa[2], filename) ~ BITS_Taxa[2],
                          grepl(BITS_Taxa[3], filename) ~ BITS_Taxa[3],
                          grepl(BITS_Taxa[4], filename) ~ BITS_Taxa[4]))
  

cols <- c(oce.colorsTemperature(30)[23],
          #oce.colorsOxygen(30)[15],
          "#E69F00",
          oce.colorsChlorophyll(20)[10],
          oce.colorsTemperature(30)[5], 
          oceColorsSalinity(20)[10], 
          oce.colorsOxygen(20)[3]
          )

relinfs_p_q1 <- relinfs %>% 
  dplyr::filter(Taxa == "Cod") %>% 
  dplyr::filter(grepl("CodAdult",filename)) %>% 
  dplyr::filter(Quarter == "Q1") %>% 
  ggplot(aes(x = rel.inf, y =  reorder_within(var, rel.inf, list(as.factor(Decade), as.factor(Quarter))), fill=var)) +
  #ggplot(aes(x = rel.inf, y = var, order = -rel.inf, fill=var)) +
  geom_bar(stat = "identity") +
  facet_grid(Decade ~ ., scale = "free")+
  theme_bw(16)+
  labs( y= "Covariate", x= "Relative Influence")+
  theme(legend.position = "none")+
  scale_y_reordered()+
  scale_fill_manual(values = rev(cols))

relinfs_p_q4 <- relinfs %>% 
  dplyr::filter(Taxa == "Cod") %>% 
  dplyr::filter(grepl("CodAdult",filename)) %>% 
  dplyr::filter(Quarter == "Q4") %>% 
  ggplot(aes(x = rel.inf, y =  reorder_within(var, rel.inf, list(as.factor(Decade), as.factor(Quarter))), fill=var)) +
  geom_bar(stat = "identity") +
  facet_grid(Decade ~ ., scale = "free_y")+
  theme_bw(16)+
  labs( y= "Covariate", x= "Relative Influence")+
  theme(legend.position = "none")+
  scale_y_reordered()+
  scale_fill_manual(values = rev(cols))

relinfs_p_q4_grid <- plot_grid(relinfs_p_q1, relinfs_p_q4, labels = c('A', 'B'), 
          label_size = 24)

ggsave("../Figures/CodAdult_relinfs_GBMs.png",relinfs_p_q4_grid, width = 24, height = 12, units = "cm", dpi = 300)

ggsave("../Figures/Cod_relinfs_nosurfoxygen_GBMs.png",relinfs, width = 18, height = 15, units = "cm", dpi = 300)



relinfs <- relinfs %>% 
  mutate(Taxasplit = case_when(grepl("CodAdult",filename) ~ "Adult Cod",
                               grepl("AdultCod",filename) ~ "Adult Cod",
                               grepl("CodJuvenile",filename) ~ "Juvenile Cod",
                               grepl("JuvenileCod",filename) ~ "Juvenile Cod",
                               TRUE ~ Taxa)) 

relinfs$Taxasplit <- factor(relinfs$Taxasplit, levels = c("Dab", "Flounder", "Plaice", "Juvenile Cod", "Adult Cod"))

relinfs %>%   
  ggplot(aes(x = var, y = (rel.inf), col = var)) +
  geom_boxplot(fill = "grey95", lwd = 0.4)+
  #geom_jitter()+
  facet_grid(Taxasplit~Quarter, scales = "free")+
  theme_bw(15)+
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 11))+
  scale_color_manual(values = rev(cols),
                     labels = c("Bottom Oxygen", 
                                "Bottom Salinity",
                                "Bottom Temperature",
                                (expression(Chlorophyll-alpha)),
                                "Mixed Layer Depth",
                                "Surface Temperature"))+
  guides(color=guide_legend(title="Variable", override.aes = list(lwd = 0.6)))+
  labs(y = "Realtive Influence (%)", x = "")
ggsave("../Figures/RelInfluence.png", width = 17, height = 15, units = "cm", dpi = 300)

  
