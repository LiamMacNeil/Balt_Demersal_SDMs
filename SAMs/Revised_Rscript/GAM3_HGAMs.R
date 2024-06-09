library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(raster)
library(viridis)

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

####################################################################################

GAM_3 <- mgcv::bam(Density_log ~ 
                    te(Latitude, Longitude, Year, by = ScientificName_WoRMS) + 
                    s(Year, ScientificName_WoRMS, bs="fs")+ 
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
                  discrete = T
)
(vRho <- acf(resid(GAM_3), plot=T))

par(mfrow=c(1,2))
acf(resid(GAM_3), main = "ACF")
pacf(resid(GAM_3), main = "pACF")
dev.off()

par(mfrow = c(2,2))
gam.check(GAM_3)
round(k.check(GAM_3),3)


#############################################################
dat$base_resid = residuals(GAM)


#Environmental clustering of residuals?
oxy <- ggplot(aes(BottomOxygen, base_resid),
              data= dat)+ 
  geom_point(aes(color=base_resid), size=0.2)+
  stat_smooth(method = "lm", linewidth=0.35)+
  geom_hline(yintercept = 0, color = "red", linetype=2)+
  facet_wrap(~Year)+
  scale_color_viridis_c(name="Residuals")+
  labs(y = "Residuals", x = expression(Oxygen~(mg~L^-1)))+
  theme_bw(6)+
  theme(legend.position = "none",
        strip.text = element_text(size=4))

sal <- ggplot(aes(BottomSalinity, base_resid),
              data= dat)+ #only look at every 2nd year
  geom_point(aes(color=base_resid), size=0.2)+
  stat_smooth(method = "lm", linewidth=0.35)+
  geom_hline(yintercept = 0, color = "red", linetype=2)+
  facet_wrap(~Year)+
  labs(x = "Salinity (\u2030)", y ="Residuals")+
  scale_color_viridis_c(name="Residuals")+
  theme_bw(6)+
  theme(legend.position = "none",
        strip.text = element_text(size=4))

temp <- ggplot(aes(bottomT, base_resid),
               data= dat)+ #only look at every 2nd year
  geom_point(aes(color=base_resid), size=0.2)+
  stat_smooth(method = "lm", linewidth=0.35)+
  geom_hline(yintercept = 0, color = "red", linetype=2)+
  facet_wrap(~Year)+
  labs(x = "  Temperature (°C)", y ="Residuals")+
  scale_color_viridis_c()+
  theme_bw(6)+
  theme(legend.position = "none",
        strip.text = element_text(size=4))

resid_env <- oxy + sal + temp + plot_layout(ncol = 1, nrow = 3)
ggsave("../Figures/Spring2024Revision/GAM_residuals_environment.png", resid_env,
       width = 17, height = 15, units = "cm" , dpi = 600)


#############################################################

summary(GAM_3)

par(mfrow = c(3,4))
plot(GAM_3, shade = T)

AIC(GAM, GAM_2, GAM_3)
BIC(GAM, GAM_2, GAM_3)

# gratia approach
one <- draw(GAM_3, residuals = F, select = 1)+
  theme_bw(12)+
  ggtitle("")

two <- draw(GAM_3, residuals = F, select = 2)+
  theme_bw(12)+
  ggtitle("")

three <- draw(GAM_3, residuals = F, select = 3)+
  #geom_line(linewidth=1.1)+
  theme_bw(12)+
  ggtitle("")

four <- draw(GAM_3, residuals = F, select = 4)+
  #geom_line(linewidth=1.1)+
  theme_bw(12)+
  ggtitle("")

five <- draw(GAM_3, residuals = F, select = 5)+
  #geom_line(linewidth=1.1)+
  theme_bw(12)+
  ggtitle("")

six <- draw(GAM_3, residuals = F, select = 6:10, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  labs(subtitle = "")+
  ggtitle("")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)

seven <- draw(GAM_3, residuals = T, select = 11, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Depth (m)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 120))+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)

eight <- draw(GAM_3, residuals = T, select = 12, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Depth (m)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 120))+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)

nine <- draw(GAM_3, residuals = T, select = 13, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Temperature (°C)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)

ten <- draw(GAM_3, residuals = T, select = 14, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Temperature (°C)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.text = element_text(size=6))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2", name = "Group") + 
  guides(color = guide_legend(override.aes = list(linewidth=2)))

eleven <- draw(GAM_3, residuals = T, select = 15, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = expression(Oxygen~(mg~L^-1)))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")

twelve <- draw(GAM_3, residuals = T, select = 16, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = expression(Oxygen~(mg~L^-1)))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")

thirteen <- draw(GAM_3, residuals = T, select = 17, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Salinity (\u2030)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        axis.title.x = element_text(size=6),
        legend.position = "none")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")

fourteen <- draw(GAM_3, residuals = T, select = 18, grouped_by = T)+
  geom_line(size=0.95)+
  theme_bw(8)+
  ggtitle("")+
  labs(x = "Salinity (\u2030)")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=6),
        plot.margin=unit(c(-0.50,0,0,0), "null"),
        legend.position = "none")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")

fifteen <- draw(GAM_3, residuals = T, select = 19, grouped_by = T)+
  theme_bw(8)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size=6),
        plot.margin=unit(c(-0.50,0,0,0), "null"))


# Env partials to collage
env <- six + seven + eight + nine +
  ten + eleven + twelve + thirteen +
  fourteen +  plot_layout(ncol = 5, nrow = 2)
ggsave("../Figures/Spring2024Revision/GAM3_env_Partials.png",env, units = "cm", width = 17, height = 15, dpi = 600)

#########
##


# Clean mapping for species-level spatiotemporal autocorrelation partials
tensor <- data.frame(one[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea_sf)

one_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(12)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~data.Year)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(two[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea_sf)

two_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(12)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~data.Year)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(three[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea_sf)

three_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(12)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~data.Year)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(four[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea_sf)

four_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(12)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~data.Year)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(five[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(icesarea_sf)

five_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(12)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  #facet_wrap(~data.Year)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA)+
  labs(x = "Longitude", y = "Latitude")

maps <- one_tens + two_tens + three_tens + four_tens + 
  five_tens +  plot_layout(ncol = 3,nrow = 2)
ggsave("../Figures/Spring2024Revision/Model_1_Partials_maps.png", maps, units = "cm", width = 17, height = 15, dpi = 600)


ggsave("../Figures/Spring2024Revision/GAM3_Dab_spatial_Partials.png",one_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/GAM3_Flounder_spatial_Partials.png",two_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/GAM3_Plaice_spatial_Partials.png",three_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/GAM3_JuvenilCod_spatial_Partials.png",four_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/GAM3_AdultCod_spatial_Partials.png",five_tens, units = "cm", width = 17, height = 15, dpi = 600)

#######################
##############################################

# Spatiotemporal prediction
##############################################
dat$pred <- predict(GAM, newdata = dat %>% 
                      dplyr::select(Latitude, Longitude,
                                    Year, ScientificName_WoRMS, 
                                    Depth,Density_log), 
                    type="response")

# Dab
map <- dat %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Limanda limanda") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.1)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.075)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf()+
  scale_color_viridis(name = expression(Predicted~log(kg~km^-2)))+
  #scale_fill_viridis(name = expression(Predicted~kg~km^-2))+
  theme_bw(9)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Dab_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Flounder
map <- dat %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.1)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.075)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf()+
  scale_color_viridis(name = expression(Predicted~log(kg~km^-2)))+
  #scale_fill_viridis(name = expression(Predicted~kg~km^-2))+
  theme_bw(9)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Flounder_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Plaice
map <- dat %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Pleuronectes platessa") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.1)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.075)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf()+
  scale_color_viridis(name = expression(Predicted~log(kg~km^-2)))+
  #scale_fill_viridis(name = expression(Predicted~kg~km^-2))+
  theme_bw(9)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Plaice_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Juvenile Cod
map <- dat %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Juvenile Gadus morhua") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.1)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.075)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf()+
  scale_color_viridis(name = expression(Predicted~log(kg~km^-2)))+
  #scale_fill_viridis(name = expression(Predicted~kg~km^-2))+
  theme_bw(9)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/JuvenileCod_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Adult Cod
map <- dat %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Adult Gadus morhua") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.1)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.075)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year)+
  coord_sf()+
  scale_color_viridis(name = expression(Predicted~log(kg~km^-2)))+
  #scale_fill_viridis(name = expression(Predicted~kg~km^-2))+
  theme_bw(9)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/AdultCod_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)
