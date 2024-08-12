rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(sf)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"

setwd(homewd)


#load the shapefile of cambodia and get the centroid
cambs <- read_sf(paste0(homewd,"/data/shapefile-climate/khm_adm1_un/khm_adm1_un.shp"))

shapefile_1 = cambs %>% st_transform(32617)
sf_cent <- st_point_on_surface(shapefile_1)

# Transform the centroids to the WGS84 CRS
sf_cent_geo <- st_transform(sf_cent, crs = 4326)

# Extract the longitude and latitude coordinates of the centroids
lon <- st_coordinates(sf_cent_geo)[,1]
lat <- st_coordinates(sf_cent_geo)[,2]

prov <- sf_cent_geo$adm1_name
lon <- st_coordinates(sf_cent_geo)[,1]
lat <- st_coordinates(sf_cent_geo)[,2]

sf_cent <- st_centroid(shapefile_1)

ggplot() + 
  geom_sf(data = shapefile_1, fill = 'white') +
  geom_sf(data = sf_cent, color = 'red') 


#store centroid lat/long 
centroid.prov <- cbind.data.frame(provname= prov, latitude=lat, longitude=lon)
centroid.prov = subset(centroid.prov, provname!="Administrative unit not available")
centroid.prov$provname[centroid.prov$provname=="Siemreap"] <- "Siem Reap"
#save 
write.csv(centroid.prov, file = paste0(homewd, "/data/centroid_provinces.csv"), row.names = F)
