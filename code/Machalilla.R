

source("code/TEAM_code.R")
machalilla.raw<-read.csv("Data/CT-PNM-2014.csv")

# then you can use ggplot2 to plot that object
library(ggplot2)
library(ggmap)
library(lubridate)
library(dplyr)
library(maptools)
# Load 'rgdal' package, which is used to read/write shapefiles and rasters
library(rgdal) 


long<-unique(machalilla.raw$longitude)
lati<-unique(machalilla.raw$latitude)
centercoord<-c(mean(subset(long, long<=1)),mean(unique(subset(lati, lati<=1))))
coordsubset<-subset(machalilla.raw,select = c(camera_trap,longitude,latitude,first_name_set_camera))

latpercam<-coordsubset %>%
  group_by(camera_trap) %>%
  summarise(latitude = median(latitude, na.rm=TRUE))

lonpercam<-coordsubset %>%
  group_by(camera_trap) %>%
  summarise(longitude = median(longitude, na.rm=TRUE))

fullcams<-dplyr::full_join(latpercam, lonpercam, by="camera_trap") #pega usando camera_trap
fullcams2<-as.data.frame(fullcams)
coordinates(fullcams2) = c("longitude", "latitude") # convierte a spatial object
# writePointsShape (fullcams2, "shp/machalilla_cams") # escribe en shapefile




map <- get_map(location = centercoord, zoom = 11, 
               source = 'google', maptype = "terrain", color ="bw")
ggmap(map) #, fullpage = TRUE)



aoi.full <- readOGR(dsn = 'shp', 'machalilla') # read shape file
aoi <- spTransform(aoi.full, CRS=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

PN_Machal<- fortify(aoi, region='name')

ap2<-ggplot(PN_Machal, aes(long, lat)) +
  geom_polygon(aes(group = group), fill = NA, colour = "green") + coord_equal()


### add camera traps
g1<-  ggmap(map, extent = 'device') + geom_point(data = machalilla.raw, 
                                                 aes(x=longitude, y=latitude), 
                                                 colour = "red", size = I(2),  
                                                 alpha = 3/8, na.rm = TRUE) 
                                                      
g2 <- ap2 + geom_point(data = machalilla.raw, 
                     aes(x=longitude, y=latitude), 
                     colour = "red", size = I(2),  
                     alpha = 3/8, na.rm = TRUE) 
  
g3<-g2 + geom_text(data = machalilla.raw, 
          aes(longitude, latitude,label=(camera_trap), size=8),
          angle = 0,
          alpha = 1/2,
          na.rm = TRUE)   
  
overlay <- stat_density2d(
  aes(x = longitude, y = latitude, fill = ..level.., alpha = ..level..),
  bins = 10, geom = "polygon",
  data = machalilla.raw) 

# con google
g3<- g1 + overlay + scale_fill_gradient(low = "blue", high = "red") + facet_wrap(~ genus, ncol = 5)
g3a<- g1 + overlay + scale_fill_gradient(low = "blue", high = "red") 

# sin google
g4<- g2 + overlay + scale_fill_gradient(low = "blue", high = "red") + facet_wrap(~ genus, ncol = 6)

############################################

#function to create binary matrices for all species at a site and sampling period. Matrix has a 1 if the species was seen in a day a 0 if not seen and NA if not sampled
#The function requires data from one sampling event and will return a list composed of 0,1 matrices, one matrix for each species.

#THIS FUNCTION WORKS WITH NEW TEAM DATA ONLY - do not use with legacy TEAM data
# this works one year at a time. Separate data in different years first


#### date fix
# unique(year(machalilla.raw$camera_trap_start_time))
machalilla.raw$photo_date2<-as.Date(machalilla.raw$photo_date, "%d-%b-%Y")
machalilla.raw$year<-year(machalilla.raw$photo_date2)
machalilla.raw$Sampling.Period<-2014
machalilla.raw$binomial<-paste(machalilla.raw$genus, machalilla.raw$specise, sep = " ")

##########################################
##### export camera points
##########################################





cam_point<-unique(machalilla.raw$camera_array)

  coordinates (machalilla.raw$longitude, machalilla.raw$latitude)

# to do
yr2011<-f.matrix.creator2(data = machalilla.raw, year = 2014) #### fix names in function





