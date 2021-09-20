### load libraries
library(leaflet)
library(tmap)
library(raster)
library(rasterVis)
library(sp)
library(terra)
library(sf)         
library(spData)        # load geographic data
library(spDataLarge)   # load larger geographic data

#library(leafletR)
library(ggplot2)
library(tidyverse)
library(tmaptools)
library(mapview)
library(mapdeck)
library(terra)
library(RColorBrewer)


## file path to ecological ecosystems
ecological_es_mean_file_path <-  "~/Desktop/AGGREGATE/Ecological_ES_Mean.tif"

## read raster files
ecological_es_mean_raster <- raster(ecological_es_mean_file_path)

### colorbrewer
display.brewer.pal(9, "Greens")
brewer.pal(9, "Greens")
raster_color <- brewer.pal(9, "YlGnBu")

pal <- colorNumeric(raster_color, domain = values(ecological_es_mean_raster) )

pal_viridis_magma <- colorNumeric("magma", values(ecological_es_mean_raster))

library(scales)
show_col(viridis_pal()(38))
