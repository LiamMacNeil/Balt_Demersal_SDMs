library(tidyverse)
library(mgcv)
library(gratia)
library(patchwork)
library(sf)
library(raster)
library(viridis)
library(itsadug)
library(oce)

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

################################################################
# Approach 1
################################################################

GAM_default <- mgcv::bam(Density_log ~ 
                   ti(Latitude, Longitude, Year_fac, d=c(2,1), bs = "fs", 
                      by=ScientificName_WoRMS, m=2) + 
                   te(Year, Depth, by = ScientificName_WoRMS, m=2, bs = "fs")+
                     s(ScientificName_WoRMS, bs="re"),
                 data = dat, 
                 family = tw(), 
                 method = 'fREML', 
                 select = T,
                 discrete=T)

rho <- acf(resid(GAM_default), 
           lag.max= 20, plot=TRUE)$acf[2]
dat <- start_event(dat, column="Year",
                   event=c("Year", "ScientificName_WoRMS"), 
                   label.event="Event")

GAM <- mgcv::bam(Density_log ~ 
                   ti(Latitude, Longitude, Year_fac, d=c(2,1), bs = "fs", 
                      by=ScientificName_WoRMS, m=2) + 
                   te(Year, Depth, by = ScientificName_WoRMS, m=2, k=c(20, 10))+
                   s(ScientificName_WoRMS, bs="re"),
                 data = dat, 
                 family = tw(), 
                 method = 'fREML', 
                 select = T,
                 discrete=T, 
                 rho = rho,
                 AR.start = start.event)
saveRDS(GAM, file = "../GAMs/GAM1.rds")
GAM <- readRDS("../GAMs/GAM1.rds")


par(mfrow = c(2,2))
gam.check(GAM)
round(k.check(GAM),3)

conc <- concrvity(GAM)
conc %>% 
  filter(.type != "worst" & .type != "observed") %>% 
  #filter(.term != "para" & .term != "s(ScientificName_WoRMS)") %>% 
  mutate(.term = gsub("ScientificName_WoRMS", "", .term)) %>% 
  #mutate(.term = gsub("te(Latitude,Longitude,Year):", "te(Lat,Lon,Year):", .term)) %>% 
  draw()+
  theme_bw(12)+
  ggtitle("")+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle=25, hjust=1, size=5))
ggsave("../Figures/Spring2024Revision/GAM_concurvity.png",
       width = 12, height = 12, units = "cm", dpi = 600)

summary(GAM)

par(mfrow = c(3,4))
plot(GAM, shade = T)

par(mfrow=c(1,2))
acf(resid(GAM_default, type = "scaled.pearson"), main = "ACF", lag.max = 20)
acf(resid(GAM, type = "scaled.pearson"), main = "ACF", lag.max = 20)
dev.off()


#####################################################################################

dat$base_resid = residuals(GAM, type = "scaled.pearson")

#Geographic clustering of residuals?
ggplot(aes(Longitude, Latitude),
       data= dat)+ #only look at every 2nd year
  geom_point(aes(color=base_resid), size=0.3, alpha=0.75)+
  #geom_sf(data=Coastline,fill="grey70", col=NA)+
  facet_wrap(~Year)+
  scale_color_viridis_c(name = "Residuals")+
  theme_bw(12)
ggsave("../Figures/Spring2024Revision/GAM_residuals_geographic.png",
       width = 17, height = 15, units = "cm" , dpi = 600)


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

################################################################################
# Partial Dependence 
################################################################################
# gratia approach
one <- draw(GAM, residuals = F, select = 1)+
  theme_bw(12)+
  ggtitle("")

two <- draw(GAM, residuals = F, select = 2)+
  theme_bw(12)+
  ggtitle("")

three <- draw(GAM, residuals = F, select = 3)+
  theme_bw(12)+
  ggtitle("")

four <- draw(GAM, residuals = F, select = 4)+
  theme_bw(12)+
  ggtitle("")

five <- draw(GAM, residuals = F, select = 5)+
  theme_bw(12)+
  ggtitle("")

six <- draw(GAM, residuals = F, select = 6, rug = F)+
  #geom_line(linewidth=1.1)+
  theme_bw(10)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"))+
  scale_y_reverse()+
  labs(caption = "",
       subtitle = expression(italic('Limanda limanda')))

seven <- draw(GAM, residuals = F, select = 7, rug = F)+
  #geom_line(linewidth=1.1)+
  theme_bw(10)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"))+
  ggtitle("")+
  scale_y_reverse()+
  labs(caption = "", 
       subtitle = expression(italic('Platichthys flesus')))

eight <- draw(GAM, residuals = F, select = 8, rug = F)+
  #geom_line(linewidth=1.1)+
  theme_bw(10)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"))+
  ggtitle("")+
  scale_y_reverse()+
  labs(caption = "",
       subtitle = expression(italic('Pleuronectes platessa')))

nine <- draw(GAM, residuals = F, select = 9, rug = F)+
  #geom_line(linewidth=1.1)+
  #scale_color_brewer(palette = "Spectral")+
  theme_bw(10)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"))+
  ggtitle("")+
  scale_y_reverse()+
  labs(caption = "",
       subtitle = expression(italic('Gadus morhua')~(Juvenile)))

#six <- draw(GAM, residuals = F, select = 6, grouped_by = T)+
#  geom_line(linewidth=1.1)+
#  scale_color_brewer(palette = "Spectral")+
#  theme_bw(12)+
#  ggtitle("")

ten <- draw(GAM, residuals = F, select = 10, rug = F)+
  theme_bw(10)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.margin=unit(c(-0.50,0,0,0), "null"))+
  ggtitle("")+
  scale_y_reverse()+
  labs(caption = "",
       subtitle = expression(italic('Gadus morhua')~(Adult)))

# Not use, indistinguishable in concurvity plots
#eleven <- draw(GAM, residuals = F, select = 11)+
#  #geom_line(linewidth=1.1)+
#  theme_bw(10)+
#  theme(
#        plot.margin=unit(c(-0.50,0,0,0), "null"))+
#  ggtitle("")+
#  labs(caption = "")
       
env <- six + seven + eight + nine + ten +  plot_layout(ncol = 3,nrow = 2)
ggsave("../Figures/Spring2024Revision/GAM_norug_env_Partials.png", env, units = "cm", width = 17, height = 15, dpi = 600)

#twelve <- draw(GAM, residuals = T, select = 12)+
  #geom_line(linewidth=1.1)+
#  theme_bw(12)+
#  ggtitle("")

# Clean mapping for species-level spatiotemporal autocorrelation partials
tensor <- data.frame(one[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(BITS_buffer)

#Center the scale on zero to avoid confusion
limit <- max(abs(tensor$data..estimate)) * c(-1, 1)

one_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(10)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~data.Year_fac)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA,
                       limit=limit)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(two[1]) %>% 
 # mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(BITS_buffer)

#Center the scale on zero to avoid confusion
limit <- max(abs(tensor$data..estimate)) * c(-1, 1)

two_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(10)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~data.Year_fac)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA, 
                       limit=limit)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(three[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(BITS_buffer)

#Center the scale on zero to avoid confusion
limit <- max(abs(tensor$data..estimate)) * c(-1, 1)

three_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(10)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~data.Year_fac)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA, 
                       limit=limit)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(four[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(BITS_buffer)

#Center the scale on zero to avoid confusion
limit <- max(abs(tensor$data..estimate)) * c(-1, 1)

four_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(10)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~data.Year_fac)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA, 
                       limit=limit)+
  labs(x = "Longitude", y = "Latitude")

tensor <- data.frame(five[1]) %>% 
  #mutate(data.Year = round(data.Year, 0)) %>% 
  mutate(Longitude = data.Longitude,
         Latitude = data.Latitude ) %>% 
  st_as_sf(coords = c("data.Longitude", "data.Latitude"),
           crs = st_crs(icesarea)) %>% 
  st_intersection(BITS_buffer)

#Center the scale on zero to avoid confusion
limit <- max(abs(tensor$data..estimate)) * c(-1, 1)

five_tens <- ggplot()+
  geom_tile(data=tensor,aes(x=Longitude, y=Latitude, fill=data..estimate))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_bw(10)+
  theme(axis.text = element_text(10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~data.Year_fac)+
  coord_sf(label_axes = "--EN",
           expand = F,
           clip = "off")+
  scale_fill_distiller(palette = "RdBu",
                       name = "Partial effect",
                       na.value = NA, 
                       limit=limit)+
  labs(x = "Longitude", y = "Latitude")

#maps <- one_tens + two_tens + three_tens + four_tens + 
#  five_tens +  plot_layout(ncol = 3,nrow = 2)
#ggsave("../Figures/Spring2024Revision/GAM_Partials_maps.png", maps, units = "cm", width = 17, height = 15, dpi = 600)

ggsave("../Figures/Spring2024Revision/Dab_spatial_Partials.png",one_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/Flounder_spatial_Partials.png",two_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/Plaice_spatial_Partials.png",three_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/JuvenilCod_spatial_Partials.png",four_tens, units = "cm", width = 17, height = 15, dpi = 600)
ggsave("../Figures/Spring2024Revision/AdultCod_spatial_Partials.png",five_tens, units = "cm", width = 17, height = 15, dpi = 600)


################################################################################
# Spatiotemporal prediction
################################################################################

dat$pred <- predict(GAM, newdata = dat %>% 
                      dplyr::select(Latitude, Longitude,Year,
                                    Year_fac, ScientificName_WoRMS, 
                                    Depth,Density_log), 
                    type="response")

# Dab
map <- dat %>% 
  mutate() %>% 
  filter(Year %in% c(2001, 2005, 2010, 2015, 2020)) %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Limanda limanda") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.25)+
  #geom_tile(aes(x = Longitude, y = Latitude, fill = pred)
  #          )+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year_fac)+
  coord_sf(expand = F)+
  scale_color_gradientn(colors = oce.colorsDensity(20),
                         name = expression(Predicted~log(kg~km^-2)),
                       na.value = NA)+
  theme_bw(10)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Dab_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Flounder
map <- dat %>% 
  filter(Year %in% c(2001, 2005, 2010, 2015, 2020)) %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Platichthys flesus") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.25)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year_fac)+
  coord_sf()+
  scale_color_gradientn(colors = oce.colorsChlorophyll(20),
                        name = expression(Predicted~log(kg~km^-2)),
                        na.value = NA)+
  theme_bw(10)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Flounder_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Plaice
map <- dat %>% 
  filter(Year %in% c(2001, 2005, 2010, 2015, 2020)) %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Pleuronectes platessa") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.25)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year_fac)+
  coord_sf()+
  scale_color_gradientn(colors = oce.colorsSalinity(20),
                        name = expression(Predicted~log(kg~km^-2)),
                        na.value = NA)+
  theme_bw(10)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/Plaice_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Juvenile Cod
map <- dat %>% 
  filter(Year %in% c(2001, 2005, 2010, 2015, 2020)) %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Juvenile Gadus morhua") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.25)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year_fac)+
  coord_sf()+
  scale_color_gradientn(colors = oce.colorsCDOM(20),
                        name = expression(Predicted~log(kg~km^-2)),
                        na.value = NA)+
  theme_bw(10)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/JuvenileCod_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

# Adult Cod
map <- dat %>% 
  filter(Year %in% c(2001, 2005, 2010, 2015, 2020)) %>% 
  st_as_sf(coords = c("x","y"), crs=st_crs(icesarea_sf)) %>% 
  filter(ScientificName_WoRMS == "Adult Gadus morhua") %>% 
  ggplot()+
  geom_sf(aes(color = pred), size=0.25)+
  #geom_tile()+
  geom_sf(data = icesarea_sf, color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.01)+
  facet_wrap(~Year_fac)+
  coord_sf()+
  scale_color_gradientn(colors = oce.colorsCDOM(20),
                        name = expression(Predicted~log(kg~km^-2)),
                        na.value = NA)+
  theme_bw(10)+
  theme(axis.text.x.bottom = element_text(angle=35, hjust=1))
ggsave("../Figures/Spring2024Revision/AdultCod_predmap.png", map,
       width = 17, height = 15, units = "cm", dpi = 600)

##############################################
