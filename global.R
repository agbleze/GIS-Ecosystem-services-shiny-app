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
library(leafsync)
library(plainview)
library(leafem)

MAPBOX = "pk.eyJ1IjoiYWdibGV6ZSIsImEiOiJja3RxYjJsdTgwNHFiMm9xZXlvazU4Z2Q3In0._u5Q5XKA-T1HCCkyzRq5iw"
## file path to ecological ecosystems
ecological_es_mean_file_path <-  "~/Desktop/AGGREGATE/Ecological_ES_Mean.tif"
#ecological_es_sum_file_path <- "~/Desktop/AGGREGATE/Ecological_ES_SUM.tif"
economic_es_mean_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ESmean_Rescale.tif"
#economic_es_mean_resam_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ES_MEAN_Clip1_Resam.tif"
social_es_mean_filepath <- "~/Desktop/AGGREGATE-social/Social_ES_Mean.tif"
economic_es_mean_clip1_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ES_MEAN_Clip1.tif"
economic_es_mean_clip1_resam_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ESmean_Rescale.tif"
#st_read("~/Desktop/AGGREGATE-social/Scocial_ES_Mean_vector.shp") -> social_es_mean_shp

tm_shape(social_es_mean_shp) + tm_polygons()
# filepath for prioritization
ecological_prioritized_filepath <- "~/Desktop/PRIORITIZE ES/Ecological_ES_Prioritized_mean.tif"
economic_prioritized_filepath <- "~/Desktop/PRIORITIZE ES/Economic_ES_Prioritized_mean.tif"
social_prioritized_filepath <- "~/Desktop/PRIORITIZE ES/SOcial_ES_Prioritized_mean.tif"
equal_weight_filepath <- "~/Desktop/PRIORITIZE ES/Equal_Weighting_All.tif"

# filepath for hotspot for prioritization
ecological_priori_hotspot_filepath <- "~/Desktop/HOTSPOT ANALYSIS/ONLY HOTSPOTS COLDSPOTS NON-SIG/Ecological_ES_Prioritized_mean_Points_Hotspot.tif"
economic_priori_hotspot_filepath <- "~/Desktop/HOTSPOT ANALYSIS/ONLY HOTSPOTS COLDSPOTS NON-SIG/Economic_ES_Prioritized_mean.tif"
social_priori_hotspot_filepath <- "~/Desktop/HOTSPOT ANALYSIS/ONLY HOTSPOTS COLDSPOTS NON-SIG/SOcial_ES_prioritized_mean_Hotspot.tif"
equal_priori_hotspot_filepath <- "~/Desktop/HOTSPOT ANALYSIS/ONLY HOTSPOTS COLDSPOTS NON-SIG/Equal_Weighting_Mean_Points_Hotspot.tif"

############### filepath to shapefiles
## file path to france shapefiles
france_shp_filepath <- "~/Desktop/FranceProjected.shp"
fran_adm3_filepath <- "~/Desktop/FRA_adm/FRA_adm3.shp"
fran_adm0_filepath <- "~/Desktop/FRA_adm/FRA_adm0.shp"
fran_adm1_filepath <- "~/Desktop/FRA_adm/FRA_adm1.shp"

## read raster files
#ecological_es_sum_raster <- raster(ecological_es_sum_file_path)
ecological_es_mean_raster <- raster(ecological_es_mean_file_path)
economic_es_mean_raster <- raster(economic_es_mean_filepath)
#economic_es_mean_resam_raster <- raster(economic_es_mean_resam_filepath)
social_es_mean_raster <- raster(social_es_mean_filepath)

economic_es_mean_clip1_raster <- raster(economic_es_mean_clip1_filepath)
economic_es_mean_clip1_resam_raster <- raster(economic_es_mean_clip1_resam_filepath)

ecological_prioritized_raster<- raster::raster(ecological_prioritized_filepath)
economic_prioritized_raster <- raster::raster(economic_prioritized_filepath)
social_prioritized_raster <- raster::raster(social_prioritized_filepath)
equal_weight_raster <- raster(equal_weight_filepath)

raster(ecological_priori_hotspot_filepath) -> ecological_priori_hotspot_raster
raster(economic_priori_hotspot_filepath) -> economic_priori_hotspot_raster
raster(social_priori_hotspot_filepath) -> social_priori_hotspot_raster
raster(equal_priori_hotspot_filepath) -> equal_priori_hotspot_raster

## read shp
france_shp <- sf::read_sf(france_shp_filepath)
france_shp_tmp <- tm_shape(france_shp) + tm_polygons()
fran_adm3_shp <- sf::read_sf(fran_adm3_filepath)
#soc_mean_hotspot_shp <- st_read(soc_mean_hotspot_shp_filepath)
fran_adm0_shp <- st_read(fran_adm0_filepath)
fran_adm1_shp <- st_read(fran_adm1_filepath)
#### create tmap objects
hotspot_economic_priori_tm <- tm_shape(economic_priori_hotspot_raster) + tm_raster()
hotspot_ecological_priori_tm <- tm_shape(ecological_priori_hotspot_raster) + tm_raster()
hotspot_social_priori_tm<- tm_shape(social_priori_hotspot_raster) + tm_raster()

## Prioritization tmap objects
economic_prioritized_tmap <- tm_shape(economic_es_mean_raster) + tm_raster()

#mapview::mapview(economic_priori_hotspot_raster)
##########################################
### colorbrewer
#display.brewer.pal(9, "Greens")
brewer.pal(9, "Greens")
raster_color <- brewer.pal(9, "YlGnBu")

pal <- colorNumeric(raster_color, domain = values(ecological_es_mean_raster) )

pal_viridis_magma <- colorNumeric("magma", values(ecological_es_mean_raster))

mapviewOptions(default = TRUE)
# mapviewGetOption("raster.palette")


m1 <- mapview(equal_priori_hotspot_raster, query.type = 'mousemove', 
              query.digits = 2, legend.opacity = 0, layer.name = "Equal prioritization") 

m2 <- mapview(economic_priori_hotspot_raster, query.type = 'mousemove',
              query.digits = 2, layer.name = "Economic prioritization")

m3 <- mapview(ecological_priori_hotspot_raster, query.type = 'mousemove',
              query.digits = 2, layer.name = "Ecological prioritization")

m4 <- mapview(social_priori_hotspot_raster, query.type = 'mousemove',
              query.digits = 2, layer.name = "Social prioritization", na.alpha = 0.1) # + paste0("Social prioritization hotspot")

sync(list(m1, m2, m3, m4)) 

#View(ecological_es_mean_raster)
#View(raster::as.data.frame(ecological_es_mean_raster))
#ecolO_na.omit <- na.omit(ecological_es_mean_raster)
#str(ecological_es_mean_raster)
#summary(ecological_es_mean_raster)
#summary(ecolO_na.omit)

# r <- raster(ncol=18,nrow=18)
# r[39:49] <- 1
# r[113:155] <- 2
# r[200] <- 6
# s <- trim(r) 
# raster::plot(r)
# raster::plot(s)
# ecolog_trim <- trim(ecological_es_mean_raster, value = NA)
# summary(ecolog_trim)
# raster::plot(ecological_es_mean_raster)
# 
# raster::minValue(social_es_mean_raster)
# raster::values(social_es_mean_raster)
# raster::Moran(social_es_mean_raster)
# soc_mean_na.omit <- na.omit(values(social_es_mean_raster$Social_ES_Mean))
#mapview(fran_adm1_shp)
