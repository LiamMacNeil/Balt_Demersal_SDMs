#::install_github("JorGarMol/VoCC", dependencies = TRUE)

library(VoCC)
library(terra)
library(oce)
library(sdmpredictors)
library(tidyverse)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

current <- stack(list.files(paste0("C:/Users/lmacneil/Documents/copernicus-processed-data/SurfaceTemp_CMEMS_2001_2020"), 
                           pattern = "\\.tif$", full.names = TRUE))

########
######## Time series plot
######## 
names(current) <- 2001:2020
plot(current[[19:20]], 
     #box = FALSE,
     axes = FALSE,
     main = paste(names(current[[19:20]])), 
     col=oce.colorsTemperature(30),
     zlim=c(-0.5,10),
     legend.args = list(text = '\u00B0C', side = 3, line = 0.5, cex=1.75),
     legend.width=1.5)

########
######## 
########

pr <- as.polygons(rast(current) > -Inf)
pr <- as(pr, "Spatial")

plot(pr)

icesarea <- gisland::read_sf_ftp("ices_areas") %>% 
  as("Spatial")
# Future layers for forecasting
predictors_fut <- list_datasets(terrestrial = F, marine = T)

predict_layers <- list_layers_future(datasets = predictors_fut)
predict_layers <- unique(predict_layers[predict_layers$dataset_code =="Bio-ORACLE",])

# Select variables trained on BART after var reduction
vars <- c("BO22_tempmean_ss")

all = c("RCP45", "RCP85")
both = c("2100")
future <- get_future_layers(vars, scenario = all, year =both)
Temp_layers <- load_layers(future$layer_code, rasterstack = T) 
Temp_layers <- crop(mask(Temp_layers, pr), extent(pr)) 
#Temp_layers <- subset(Temp_layers, order(c(1,3,2,4)))

currentresampled <- projectRaster(current, Temp_layers, method = 'bilinear')

Temp_layers_all <- stack(currentresampled, Temp_layers)

temp <- stack(currentresampled, Temp_layers[[2]])
# monthly to annual averages
vt_balt <- tempTrend(currentresampled, th = 20)
#vt_balt <- tempTrend(temp, th = 10)

pdf("Surface_Veloicty.pdf")
plot(vt_balt[[1]], box=F, axes=F, col=oce.colorsTemperature(20))
dev.off()
# spatial gradient
vg <- spatGrad(currentresampled, projected = FALSE, th=0.005)

# gradient-based climate velocity
gv <- gVoCC(vt_balt, vg)
gv <- crop(mask(gv, pr), extent(pr))
crs(gv) <- "+proj=longlat +datum=WGS84"
writeRaster(gv[[1]], filename = paste0(names(gv)[1], ".tif"), overwrite=T)

gv$voccMag[] <- ifelse(is.infinite(gv$voccMag[]), NA, gv$voccMag[]) # replace inf with NAs (just in case)
vel <- gv[[1]]
ang <- gv[[2]]

mn <- calc(Temp_layers_all, mean, na.rm = T)

x1 <- crop(gv[[1]], extent(gv))
x2 <- crop(gv[[1]], extent(gv))   
extent(x1) <- extent(gv)
velc <- merge(x1, x2)
# crop to the desired extent 

lonlat <- data.frame(xyFromCell(velc, 1:ncell(velc)))
lonlat$vel <- raster::extract(vel, lonlat)
lonlat$ang <- raster::extract(ang, lonlat[,1:2])
lonlat$mn <- raster::extract(mn, lonlat[,1:2])
lonlat <- na.omit(lonlat)

mnd <- disaggregate(mn, 4)
veld <- disaggregate(vel, 4)
angd <- disaggregate(ang, 4)
lonlat <- na.omit(data.frame(xyFromCell(veld, 1:ncell(veld)), veld[], angd[], mnd[]))[,1:2]

lonlat <- na.omit(data.frame(xyFromCell(vel, 1:ncell(vel)), vel[], ang[], mn[]))[,1:2]
traj <- voccTraj(lonlat, vel, ang, mn, tyr = 80, correct = TRUE)

lns <- trajLine(x = traj, projx = "+proj=longlat +datum=WGS84")

pdf("veltrajlines.pdf")
plot(pr, axes=F)
plot(lns, add = TRUE, lwd=0.01, col=oce.colorsTemperature(25))
dev.off()

pdf("VelMag.pdf")
plot(gv[[1]], box=F, axes=F, col=oce.colorsTemperature(25),
     legend.args = list(text = expression(km~yr^-1), side = 3, line = 0.5, cex=1.75),
     legend.width=1.5)
#plot(lns, add = TRUE, lwd=1e-10)
dev.off()

# Raw trajClass funtion from VoCC (need to adjust rasterstack "b")

trajClas <- function(traj, vel, ang, mn, trajSt, tyr, nmL, smL , Nend, Nst, NFT, DateLine = FALSE){
  
  TrajEnd <- TrajFT <- TrajSt <- IntS <- BounS <- TrajClas <- raster(ang)
  
  # add cell ID to the data frame
  traj <- data.table(traj)
  traj$cid <- cellFromXY(ang, traj[,1:2])
  
  # A. Number of traj starting from each cell
  TrajSt[!is.na(ang[])] <- 10
  
  # B. Number of traj ending in each cell
  tr <- traj[, .SD[.N], by = trajIDs]  # subset last point of each trajectory
  enTraj <- tr[,.N, by = cid]
  TrajEnd[!is.na(vel)] <- 0
  TrajEnd[enTraj$cid] <- enTraj$N
  
  # C. Number of traj flowing through each cell
  cxtrj <- unique(traj, by = c("trajIDs", "cid"))
  TotTraj <- cxtrj[,.N, by = cid]       # total number of trajectories per cell
  TrajFT[!is.na(vel)] <- 0
  TrajFT[TotTraj$cid] <- TotTraj$N
  TrajFT <- TrajFT - TrajEnd - TrajSt   # subtract traj starting and ending to get those actually transversing the cell
  TrajFT[TrajFT[] < 0] <- 0   # to avoid negative values in ice covered cells
  
  # C. Identify cell location for internal sinks (groups of 4 endorheic cells with angles pointing inwards)
  ll <- data.table(xyFromCell(ang, 1:ncell(ang)))
  ll[,1:2] <- ll[,1:2] + 0.1   # add small offset to the centroid
  a <- fourCellsFromXY(ang, as.matrix(ll[,1:2]))
  #########################  This line is throwing error for b raster below (dataframe not list needed)
  #a <- t(apply(a, 1, sort))
  # If date line crossing, correct sequences on date line
  if(DateLine == TRUE){
    a[seq(ncol(ang), by = ncol(ang), length = nrow(ang)),] <- t(apply(a[seq(ncol(ang), by = ncol(ang), length = nrow(ang)),], 1, function(x) {x[c(2,1,4,3)]}))
  }
  # Extract the angles for each group of 4 cells
  b <- matrix(raster::extract(ang, as.vector(a)), nrow = ncell(ang), ncol = 4, byrow = FALSE)
  ll[, c("c1", "c2", "c3", "c4", "ang1", "ang2", "ang3", "ang4") := data.frame(a, b)]
  # now look if the 4 angles point inwards (internal sink)
  ll[, c("d1", "d2", "d3", "d4") := .(((ang1 - 180)  *  (90 - ang1)), ((ang2 - 270)  *  (180 - ang2)), ((ang3 - 90)  *  (0 - ang3)), ((ang4 - 360)  *  (270 - ang4)))]
  ll[, isink := 0L]
  ll[d1 > 0 & d2 > 0 & d3 > 0 & d4 > 0, isink := 1L]
  # get the cids for the cells contained in the sinks
  celdas <- ll[isink == 1, 3:6]
  IntS[!is.na(vel)] <- 0
  IntS[c(celdas[[1]], celdas[[2]], celdas[[3]], celdas[[4]])] <- 1
  
  # D. Identify cell location for boundary sinks (coastal cells which are disconected from cooler climates under warming or warmer climates under cooling)
  # detect coastal cells
  coast <- suppressWarnings(boundaries(vel, type = 'inner', asNA = TRUE))       # to avoid warning for coercing NAs via asNA = TRUE
  # make a list of vel values and SST values for each coastal cells and their marine neighbours
  cc <- na.omit(data.table(cid = 1:ncell(vel), coast = coast[]))
  ad <- adjacent(vel, cc$cid, 8, sorted = TRUE, include = TRUE) # matrix with adjacent cells
  ad <- data.table(coastal = ad[,1], adjacent = ad[,2], cvel = vel[ad[,1]], ctemp = mn[ad[,1]], atemp = mn[ad[,2]])
  # locate the sinks
  ad <- na.omit(ad[ad$cvel != 0,])      # remove cells with 0 velocity (ice) and with NA (land neighbours)
  j <- ad[, ifelse(.SD$cvel > 0, all(.SD$ctemp <= .SD$atemp), all(.SD$ctemp >= .SD$atemp)), by = coastal]
  setkey(j)
  j <- unique(j)
  BounS[!is.na(vel)] <- 0
  BounS[unique(subset(j$coastal, j$V == TRUE))] <- 1
  
  # Total number of trajectories per cell and proportions per cell
  TrajTotal <- calc(stack(TrajSt, TrajFT, TrajEnd), sum, na.rm = TRUE)
  TrajTotal[is.na(ang[])] <- NA
  PropTrajEnd <- (TrajEnd/TrajTotal)*100
  PropTrajFT <- (TrajFT/TrajTotal)*100
  PropTrajSt <- (TrajSt/TrajTotal)*100
  
  # reclassify by traj length
  rclM <- matrix(c(0, (nmL/tyr), 1, (nmL/tyr), (smL/tyr), 2, (smL/tyr), Inf, 3), ncol=3, byrow=TRUE)
  v <- raster(vel)
  v[] <- abs(vel[])
  ClassMov <- reclassify(v, rclM)
  
  # Classify the cells
  TrajClas[!is.na(vel)] <- 0
  # capture non-moving (1)
  TrajClas[ClassMov[] == 1] <- 1
  # capture slow-moving (2)
  TrajClas[ClassMov[] == 2] <- 2
  # capture internal (3) and (4) boundary sinks
  TrajClas[IntS[] == 1] <- 3
  TrajClas[BounS[] == 1] <- 4
  # Classify remaining cells into sources(5), rel sinks(6), corridors(7), divergence(8) and convergence(9)
  d <- na.omit(data.table(cid = 1:ncell(TrajClas), val = TrajClas[]))
  d <- d[val == 0, 1]
  d[, Nend := PropTrajEnd[d$cid]]
  d[, Nst := PropTrajSt[d$cid]]
  d[, NFT := PropTrajFT[d$cid]]
  d$clas <- ifelse(d$Nend == 0, 5, ifelse(d$Nend > Nend & d$Nst < Nst, 6, ifelse(d$NFT > NFT, 7,ifelse(d$Nend < d$Nst, 8, 9))))
  TrajClas[d$cid] <- d$clas
  
  # return raster
  s <- stack(PropTrajEnd, PropTrajFT, PropTrajSt, ClassMov, IntS, BounS, TrajClas)
  names(s) <- c("PropEnd", "PropFT", "PropSt", "ClassL", "IntS", "BounS", "TrajClas")
  return(s)
}

clas <- trajClas(traj, vel, ang, mn, trajSt = 10, tyr = 80, nmL = 50, smL = 150,
                 Nend = 45, Nst = 15, NFT = 70, DateLine = FALSE)

# Define first the colour palette for the full set of categories
my_col = c('forestgreen', 'darkseagreen1', 'maroon')

#movement pace
my_col <- my_col[sort(unique(clas[[4]][]))]

r <- ratify(clas[[4]])
rat_r <-levels(r)[[1]]
rat_r$trajcat <- c("N-M", "S-M", "F-M")[sort(unique(clas[[4]][]))]
levels(r) <- rat_r
# Produce the plot using the rasterVis levelplot function
rasterVis::levelplot(r, col.regions = my_col, xlab = NULL, ylab = NULL, scales = list(draw=FALSE))

# Keep only the categories present in our raster
my_col = c('grey50', 'forestgreen', "black",'darkorange1', 'mediumblue', 'magenta3',
           'purple', 'cyan3', 'yellow')

my_col <- my_col[sort(unique(clas[[7]][]))]

# Classify raster / build attribute table
r <- ratify(clas[[7]])
rat_r <-levels(r)[[1]]
rat_r$trajcat <- c("Non-Moving", "Slow-Moving", "Internal Sink",
                   "Boundary Sink", "Source", "Relative Sink", "Corridor", 
                   "Divergence", "Convergence")[sort(unique(clas[[7]][]))]
levels(r) <- rat_r
# Produce the plot using the rasterVis levelplot function
rasterVis::levelplot(r, col.regions = my_col, xlab = NULL, ylab = NULL, scales = list(draw=FALSE),
                     colorkey=list(labels=list(cex=1.5)))

#(N-M) non-moving, short trajectories, little thermal shift, migration unexpected, stable
#(S-M) slow-moving, moderate trajectory length, slow/moderate thermal shift, little migration, stable
#(IS) internal sinks, 
#(BS) boundary sinks, 
#(Srce) sources, 
#(RS) relative sinks, 
#(Cor) corridors, 
#(Div) divergence and 
#(Con) convergence.

for (p in 1:nlayers(clas)) {
  writeRaster(clas[[p]], filename = paste0(names(clas)[p],"_Velocity", ".tif"))
}


################ Distance-based approach currently unused

# start with distance-based approach because we have averages representing periods

delta <- stack(mean(Temp_layers_all[[1:20]]), na.rm = TRUE,Temp_layers_all[[21]], na.rm = TRUE ,Temp_layers_all[[22]], na.rm = TRUE)

pdf("SurfaceTemp_Comapre.pdf")
plot(delta, box=F, axes=F, col=oce.colorsTemperature(50), zlim=c(0,15), main="")
dev.off()

# prepare the data frame with the necessary variables
clim <- na.omit(data.frame(getValues(delta), cid = 1:ncell(Temp_layers_all)))
clim[,c("x","y")] <- xyFromCell(Temp_layers_all, clim$cid)

dist_v_balt <- dVoCC(clim, n = 1, tdiff = 80, method = "Single", climTol =  10, 
           geoTol = 100, distfun = "GreatCircle", trans = NA, lonlat = T)

# Change sign as needed and create the distance-based velocity raster
#v_balt <- na.omit(v_balt)
ind <- which(delta[[1]][dist_v_balt$focal] > delta[[2]][dist_v_balt$target])
dist_v_balt$velBis <- dist_v_balt$vel
dist_v_balt$velBis[ind] <- dist_v_balt$vel[ind] * -1
# put output in raster format
dv <- raster(gv)
dv[dist_v_balt$focal] <- dist_v_balt$velBis

'''
BO_2050_RCP45_layers <- raster("C:/Users/lmacneil/Documents/outputs/BO_Cod/BO_2050_RCP45/BO22_RCP45_2050_tempmean_ss.tif")

BO_2050_RCP85_layers <- raster("C:/Users/lmacneil/Documents/outputs/BO_Cod/BO_2050_RCP85/BO22_RCP85_2050_tempmean_ss.tif")

BO_2100_RCP45_layers <- raster("C:/Users/lmacneil/Documents/outputs/BO_Cod/BO_2100_RCP45/BO22_RCP45_2100_tempmean_ss.tif")

BO_2100_RCP85_layers <- raster("C:/Users/lmacneil/Documents/outputs/BO_Cod/BO_2100_RCP85/BO22_RCP85_2100_tempmean_ss.tif")
'''
