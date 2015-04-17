


# load library
library(ggmap)
library(rgdal)

pacoche<-read.csv("Data/Pacoche.csv" ,header = TRUE, fileEncoding ="latin3")
pacoche<-pacoche[-28,] # delete one NA

pacoche2<-read.csv("Data/Pacoche2.csv" ,header = TRUE, fileEncoding ="latin3")

machalilla<-read.csv("Data/Machalilla.csv" ,header = TRUE, fileEncoding ="latin3")


# calculate midle point
Mach_meanLat<-mean(machalilla$Lat, na.rm = T) 
Mach_meanLong<-mean(machalilla$Lon, na.rm = T)
Mach_location<-c(Mach_meanLong, Mach_meanLat)

pacoche2_meanLat<-mean(pacoche2$Lat, na.rm = T) 
pacoche2_meanLong<-mean(pacoche2$Lon, na.rm = T)
pacoche2_location<-c(pacoche2_meanLong, pacoche2_meanLat)


#Plot all
plot(x=pacoche$Lon,y=pacoche$Lat,asp=1) #Boring !
plot(x=machalilla$Lon,y=machalilla$Lat,asp=1) #Boring !

## more elaborated map using ggmap
require(ggmap) # same as ggplot2 but for maps

###################################
#### Machalilla ······
##################################

wmap.base.m = qmap(Mach_location, zoom = 12) #Get a base map from google with location in the center

wmap.base.m 

wmap.m<- wmap.base.m + geom_point(data = machalilla, 
                           aes(Lon, Lat), 
                           size = I(3), 
                           color = "red", na.rm = TRUE)
#Plot plus text
wmap.m + geom_text(data = machalilla, 
                  aes(Lon, Lat,label=(Code), size=8),
                  angle = 0) 


###################################
#### Pacoche ······
##################################

wmap.base.p = qmap(pacoche2_location, zoom = 12) #Get a base map from google with location in the center

wmap.base.p

wmap.p<- wmap.base.p + geom_point(data = pacoche2, 
                           aes(Lon, Lat), 
                           size = I(3), 
                           color = "red", na.rm = TRUE)
#Plot plus text
wmap.p + geom_text(data = pacoche2, 
                  aes(Lon, Lat,label=(Code), size=8),
                  angle = 0,
                  alpha = 1/2,
                  na.rm = TRUE) 
##########################

replace<-subset(pacoche2, NOVALIDA=="X")
index<-which(is.na(replace$Lat))
replace2<-replace[-index,]

wmap.base.p2 = qmap(pacoche2_location, zoom = 13, maptype = 'hybrid') #Get a base map from google with location in the center

wmap.p2<- wmap.base.p2 + geom_point(data = replace2, 
                                  aes(Lon, Lat), 
                                  size = I(3), 
                                  color = "red", na.rm = TRUE)
wmap.p2 + geom_text(data = replace2, 
                   aes(Lon, Lat,label=(Code),  size=8),
                   angle = 0,
                   alpha = 1/2,
                   na.rm = TRUE) 

library(maptools)
library(rgdal)
coordinates(replace2)<- c("Lon", "Lat")
## In order for the points to be displayed in the correct place they need to be re-projected to WGS84 geographical coordinates.
geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(replace2)<-geo
replace2_wgs84<- spTransform(replace2, CRS= geo)
writeOGR(replace2_wgs84, dsn="replace_pacoche.kml", layer= "replace2_wgs84", driver="KML", dataset_options=c("NameField=code"))



#########################
## Changing projection ##
#########################
require(rgdal)

### take a look to http://spatialreference.org
latlon<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

campoints.m1<-as.data.frame(cbind(machalilla$Lon,machalilla$Lat)) #extract the points from data
campoints.m<-cbind(campoints.m1, machalilla$Code) 

campoints.p1<-as.data.frame(cbind(pacoche$Lon,pacoche$Lat)) #extract the points from data
campoints.p<-cbind(campoints.p1, pacoche$Code) 

colnames(campoints.m)<-c("lon","lat", "cam")
colnames(campoints.p)<-c("lon","lat", "cam")
# change from data frame to Object of class SpatialPoints
coordinates(campoints.m)<-c("lon","lat") # 
proj4string(campoints.m)<- latlon #Make spatial points

coordinates(campoints.p)<-c("lon","lat") # 
proj4string(campoints.p)<- latlon #Make spatial points

library(maptools)
writePointsShape(campoints.p, "shp/campoints_p") # export to shp
writePointsShape(campoints.m, "shp/campoints_m")



