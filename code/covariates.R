

# then you can use ggplot2 to plot that object
library(ggplot2)
library(ggmap)
library(lubridate)
library(dplyr)
library(maptools)
# Load 'rgdal' package, which is used to read/write shapefiles and rasters
library(rgdal) 
library(dismo)
library(sdm)

machalilla.raw<-read.csv("Data/CT-PNM-2014.csv") # ojo falta la camara 3-12

long<-unique(machalilla.raw$longitude)
lati<-unique(machalilla.raw$latitude)
centercoord<-c(mean(subset(long, long<=1)),mean(unique(subset(lati, lati<=1))))
coordsubset<-subset(machalilla.raw,select = c(camera_trap,longitude,latitude,first_name_set_camera))

# get elevation
elevation<-getData(name = "SRTM",lon=centercoord[1], lat=centercoord[2])


cam.cords<-as.data.frame(distinct(coordsubset))
coordinates(cam.cords) <- ~longitude+latitude #make sppatial data fram
geo <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #def cord
proj4string(cam.cords)<-geo # set cords

e<-extent (-80.9, -80.5, -1.75,-1.35)
elevation.crop<-crop(elevation, e)

plot(elevation.crop)
plot(cam.cords, add=T)
slope<-terrain(elevation.crop, opt='slope', unit='degrees', neighbors=4)


cam.cords.sp<-SpatialPoints(cam.cords)
proj4string(cam.cords.sp)<-geo 
# etract values
elev.ovr <- extract(elevation.crop, cam.cords, method='bilinear')
slope.ovr <- extract(slope, cam.cords, method='bilinear')

# add to table
cam.cords$elev<-elev.ovr
cam.cords$slope<-slope.ovr




