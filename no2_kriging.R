options(scipen = 100) # To see long decimal points
memory.size() # for WindowsOS
memory.limit(99999) # for WindowsOS

library(tidyverse)
library(sf)
library(raster)
library(rgdal)
library(automap)
library(gridExtra)


# set working directory so I know where the .zip file will be located
getwd()
#setwd(dir = "/some/path/")

# on the GitHub repository of interest
download.file(url = "https://github.com/mrsensible/GISRUK2019/archive/master.zip", 
              destfile = "GISRUK2019-master.zip")

# unzip the .zip file
unzip(zipfile = "GISRUK2019-master.zip")

# examine the contents
list.files('./GISRUK2019-master')
list.files('./GISRUK2019-master/data')

# Set Workding Directory
setwd('./GISRUK2019-master')


# Load NO2 Pollution data
load("data/no2_jan.RData")

# Import moritoring stations from Seoul
stations <- read_sf("data/stations_10km.shp")
stations_df <- stations %>% st_set_geometry(NULL)

# Import Seoul Shapefile
seoul <- read_sf("data/Seoul_City.shp") %>% as('Spatial') %>% fortify()

no2.winter <- merge(no2.win.12hr, stations_df, by.x = c("Station.ID", "X", "Y"), by.y = c("Station", "X", "Y"))
coordinates(no2.winter) <- ~X+Y
proj4string(no2.winter) <- CRS("+init=epsg:5181")


#--Put Multiple Plots on a Single Page in R with spplot--##
plots <- lapply(names(no2.winter)[3:22], function(.x) spplot(no2.winter,.x))
do.call(grid.arrange,plots)


#################################################################
#--Generate auto Semivariograms in need to create Kriging maps--#
#################################################################
myList <- list()

for(i in 1:20) { 
  myList[[length(myList)+1]] <- autofitVariogram(no2.winter[[i+2]] ~ 1, no2.winter)
}
semvar <- lapply(myList, function(x) plot(x))
do.call(grid.arrange, semvar[1:4])


### Create gridcells for interpolation
seoul_grid <- data.frame(expand.grid(X = seq(min(no2.winter$X), max(no2.winter$X), length=200),
                                     Y = seq(min(no2.winter$Y), max(no2.winter$Y), length=200)))
coordinates(seoul_grid) <- ~X+Y
proj4string(seoul_grid) <- CRS("+init=epsg:5181") #Korean Central Belt 2000


##############
#--Kriging--##
##############
sum.squares <- vector()
var.model <- data.frame()
pred.model <- seoul_grid@coords


# This iteration takes 5 minutes!!

for(i in 1:20) {
  kriging_new <- autoKrige(no2.winter@data[,i+2]~ X + Y,
                           nmax = 20000,
                           input_data = no2.winter, 
                           new_data = seoul_grid)
  sum.squares <- append(sum.squares, kriging_new$sserr)
  kriging_new$var_model <- data.frame(y=i,kriging_new$var_model)
  var.model <- rbind(var.model, kriging_new$var_model)
  xyz <- as.data.frame(kriging_new$krige_output)
  p <- data.frame(xyz[,'var1.pred'])
  colnames(p) <- colnames(no2.winter@data)[i+2]
  pred.model <- cbind(pred.model, p)
} 

##-- Add ColNames
colnames(pred.model) <- c("X", "Y", "jan01d", "jan01n", "jan02d", "jan02n","jan03d", "jan03n", "jan04d", "jan04n", "jan05d", "jan05n", "jan06d", "jan06n", "jan07d", "jan07n", "jan08d", "jan08n", "jan09d", "jan09n", "jan10d", "jan10n")


##-- Mean and variance to display on map
stat <- pred.model %>% dplyr::select(-c(X,Y)) %>% 
        gather(factor_key = T) %>% 
        group_by(key) %>% summarise(mean= round(mean(value),1), sd= round(sd(value),1), 
                                    max = max(value),min = min(value)) %>% 
        rename(Hour = key)

##############################################
##-- Final Map: Kriging Interpolation map --##
##############################################

ras.krige.df <- pred.model %>% 
  reshape2::melt(id = c("X", "Y"), variable.name = "Hour", value.name = "NO2") 

ras.krige.df %>% 
  ggplot() +
  geom_tile(aes(x = X, y = Y, fill = NO2)) +
  scale_fill_distiller(palette = "Spectral", na.value = NA, limits = c(0,125), breaks = c(0,25,50,75,100,125)) +
  geom_contour(aes(x = X, y = Y, z = NO2),bins = 20, colour = "grey40", alpha = 0.7) +
  scale_color_gradientn(limits = c(100,300), colours=c("orangered", "firebrick")) +
  geom_path(data = seoul, aes(x = long, y = lat), color = 'black', size = 1) +
  geom_text(data = stat, aes(x = 187000,  y = 434000, label = paste0("mean = " , mean)), size = 3) + 
  geom_text(data = stat, aes(x = 184000,  y = 430500, label = paste0("sd = " , sd)), size = 3) + 
  labs(title = "Kriging Interpolation for NO2 Mapping: An example of Seoul", 
       subtitle = "Hourly data aggregated to Days and Nights") +
  facet_wrap(~ Hour, ncol = 8) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 20),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15)                                  
  ) -> final # 1200 x 550 

# Export PNG
png("plot.png", width=1200, height=550, res=100)
final
dev.off()


# RMSE
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))}

for (i in 3:length(pred.model)){
  RMSE(mean(pred.model[, i]), pred.model[, i]) %>% print()
}


# convert to Raster Bricks
krige <- rasterFromXYZ(pred.model, 
                       crs="+proj=tmerc +lat_0=38 +lon_0=127 +k=1 +x_0=200000 +y_0=500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                       digits=5)

# Write Raster
writeRaster(krige, filename="seoul_no2_multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


