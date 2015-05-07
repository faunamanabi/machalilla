




###############################################################################
# This script shows an example of how to use the "gfcanalysis" R package put 
# together by Alex Zvoleff (azvoleff@conservation.org) for working with the 
# Hansen et al. 2013 Global Forest Change dataset. Contact Alex if you notice 
# any issues or have problems using the package.
#
# See the help files for the functions below for more information. For example, 
# type "?download_tiles" in R to see the help file for the "download_tiles" 
# function.
# 
# NOTE: the gfcanalysis package must be installed before this script will run.  
# Run the "install_gfcanalysis.R" script to install/update the gfcanalysis 
# package.
###############################################################################

# Load the gfcanalysis package
library(gfcanalysis)
# Load 'rgdal' package, which is used to read/write shapefiles and rasters
library(rgdal)

# Start a cluster for parallel processing. You must install the "snow" package 
# for this to work. Comment out this line, and the endCluster() line at the end 
# of this script, if you do NOT want gfcanalysis to run in parallel.
if (require(snow)) beginCluster()

# Indicate where we want to save GFC tiles downloaded from Google. For any 
# given AOI, the script will first check to see if these tiles are available 
# locally (in the below folder) before downloading them from the server - so I 
# recommend storing ALL of your GFC tiles in the same folder. For this example 
# we will save files in the current working directory folder.
output_folder <- 'HansenData'

# Set the threshold for forest/non-forest based on the treecover2000 layer in 
# the GFC product
forest_threshold <- 75

###############################################################################
# Download data from Google server for a given AOI
###############################################################################

# Load a demo AOI from the P drive - notice that first we specify the folder 
# the shapefile is in, and then the name of the shapefile without the '.shp'
# aoi <- readOGR(system.file('extdata', package='gfcanalysis'), 'ZOI_NAK_2012')
# aoi <- readOGR('.', 'ZOI_NAK_2012_EEsimple')) # . is current directory

aoi.full <- readOGR(dsn = 'shp', 'pacoche')
# aoi.full2 <- spTransform(aoi.full, CRS=CRS("+proj=utm +zone=18 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
aoi <- spTransform(aoi.full, CRS=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# aoi <- aoi.full2[aoi.full2$Unique_ID == 26,] # first square


# Calculate the google server URLs for the tiles needed to cover the AOI
tiles <- calc_gfc_tiles(aoi)
print(length(tiles))

# just to chk tile vs map
plot(tiles)
plot(aoi, add=TRUE, lty=2, col="#00ff0050") # full map
# plot(aoi.full.geo[aoi.full.geo$Unique_ID == 26,], add=TRUE, lty=2, col="red") # tile 1

# Check to see if these tiles are already present locally, and download them if 
# they are not.
download_tiles(tiles, output_folder)

# Extract the GFC data for this AOI from the downloaded GFC tiles, mosaicing 
# multiple tiles as necessary (if needed to cover the AOI).
gfc_data <- extract_gfc(aoi, output_folder)

# Save the output data to a GeoTIFF (can also save in ENVI format, Erdas 
# format, etc.)
gfc_data <- writeRaster(gfc_data, filename='data\\pacoche.tif', overwrite=TRUE)

###############################################################################
# Performing thresholding and calculate basic statistics
###############################################################################

# Calculate and save a thresholded version of the GFC product
gfc_thresholded <- threshold_gfc(gfc_data, forest_threshold=forest_threshold, 
                                 filename="data\\pacoche_extract_thresholded.tif", overwrite=TRUE)
# see http://azvoleff.com/articles/analyzing-forest-change-with-gfcanalysis/
# for threshold details

# change gfc_thresholded to UTM to get meters
projUTM <- CRS('+proj=utm +zone=16')
gfc_thresholded_utm<-projectRaster(gfc_thresholded, crs=projUTM) 
# Calculate annual statistics on forest loss/gain
gfc_stats <- gfc_stats(aoi, gfc_thresholded_utm)

# Save statistics to CSV files for use in Excel, etc.
write.csv(gfc_stats$loss_table, file='data\\pacoche_losstable.csv', row.names=FALSE)
write.csv(gfc_stats$gain_table, file='data\\pacoche_gaintable.csv', row.names=FALSE)

###############################################################################
# Make visualization of forest change
###############################################################################

# Calculate and save a thresholded annual layer stack from the GFC product 
# (useful for simple visualizations, etc.)
gfc_thresholded_annual <- annual_stack(gfc_thresholded)
writeRaster(gfc_thresholded_annual, filename='data\\pacoche_thresholded_annual.tif', overwrite=TRUE)

# Save a simple visualization of the thresholded annual layer stack (this is 
# just an example, and is using the data in WGS84. The data should be projected 
# for this).
animate_annual(aoi, gfc_thresholded_annual)

# Stop the parallel processing cluster
if (require(snow)) endCluster()

# my.at <- seq(0, 6, 1)
# 
# myColorkey <- list(at=my.at, ## where the colors change
#                    labels=list(
#                      at=my.at ## where to print labels
#                    ))
# levelplot(r, at=my.at, colorkey=myColorkey)

source("multiplot.R")
library(ggplot2)
size_scale = 1
maxpixels = 50000

library(plyr)
pyr<-list()
aoi_tr <- spTransform(aoi.full, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
aoi_points <- fortify(aoi_tr, region = "id")
aoi_df <- join(aoi_points, aoi@data, by = "id")

# paco<- ggplot(aoi_df) + 
#   aes(long,lat,group=group) + 
#   geom_polygon(alpha = 0.1) +
#   geom_path(linetype = label, alpha = 0.7) +
#   coord_equal() 

for (i in 1:12){
  p<-gplot(gfc_thresholded_annual[[i]], maxpixels = 50000,extent = 'device', legend = 'topleft') + 
    geom_tile(aes(fill = factor(value, levels = c(1, 2, 3, 4, 5, 6, 0)))) + coord_fixed() + 
    scale_fill_manual("Cover", breaks = c("1", "2", "3", 
                                          "4", "5", "6", "0"), # labels = c("Forest", "Non-forest", "Forest loss", 
                                                                         # "Forest gain", "Loss and gain", "Water", 
                                                                          #"No data"), 
                                                              values = c("#008000", "#ffffb3", "#ff0000", 
                                                                          "#0000ff", "#ff00ff", "#c0c0c0", 
                                                                          "#101010")) +
    theme(legend.position = "none")
  
  
  
  mp<-p + geom_path(data = aoi_df, aes(x = long, y = lat, group = group, alpha = 0.05), ) + 
    labs(
    #colour = "Pacoche",
    title = names(gfc_thresholded_annual[[i]]), 
    x = "Lon", y= "Lat") +
    scale_alpha(range = c(.4, .75), guide = FALSE)
    
  
  pyr[[i]] <- mp
}


########################## end loop
p<-gplot(gfc_thresholded_annual[[13]], maxpixels = 50000,extent = 'device', legend = 'topleft') + 
  geom_tile(aes(fill = factor(value, levels = c(1, 2, 3, 4, 5, 6, 0)))) + coord_fixed() + 
  scale_fill_manual("Cover", breaks = c("1", "2", "3", 
                                        "4", "5", "6", "0"),  labels = c("Forest", "Non-forest", "Forest loss", 
                     "Forest gain", "Loss and gain", "Water", 
                    "No data"), 
                    values = c("#008000", "#ffffb3", "#ff0000", 
                               "#0000ff", "#ff00ff", "#c0c0c0", 
                               "#101010")) 
  


mp<-p + geom_path(data = aoi_df, aes(x = long, y = lat, group = group, alpha = 0.05), ) + 
  labs(
    #colour = "Pacoche",
    title = names(gfc_thresholded_annual[[13]]), 
    x = "Lon", y= "Lat") +
  scale_alpha(range = c(.4, .75), guide = FALSE) 


pyr[[13]] <- mp

#####################################


multiplot(pyr[[1]], pyr[[2]], pyr[[3]], pyr[[4]],pyr[[5]], pyr[[6]],
          pyr[[7]], pyr[[8]], pyr[[9]], pyr[[10]],pyr[[11]], pyr[[12]], pyr[[13]],
                    cols=3)

multiplot(plotlist = pyr, cols=4)
########################################################

library(gtable)

## set up individual plots

# map
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$legend.position <- "none"
new_theme_empty$plot.margin <- structure(c(0, 0, 0, 0), unit = "lines", valid.unit = 3L, class = "unit")
new_theme_empty$axis.ticks <- element_blank()
new_theme_empty$axis.title.x <- element_blank()
new_theme_empty$axis.title.y <- element_blank()

















