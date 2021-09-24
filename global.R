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
library(rasterDT)

MAPBOX = "pk.eyJ1IjoiYWdibGV6ZSIsImEiOiJja3RxYjJsdTgwNHFiMm9xZXlvazU4Z2Q3In0._u5Q5XKA-T1HCCkyzRq5iw"
## file path to ecological ecosystems
ecological_es_mean_file_path <-  "~/Desktop/AGGREGATE/Ecological_ES_Mean.tif"
#ecological_es_sum_file_path <- "~/Desktop/AGGREGATE/Ecological_ES_SUM.tif"
economic_es_mean_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ESmean_Rescale.tif"
#economic_es_mean_resam_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ES_MEAN_Clip1_Resam.tif"
social_es_mean_filepath <- "~/Desktop/AGGREGATE-social/Social_ES_Mean.tif"
economic_es_mean_clip1_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ES_MEAN_Clip1.tif"
# economic_es_mean_clip1_resam_filepath <- "~/Desktop/AGGREGATE ECONOMIC/Economic_ESmean_Rescale.tif"
#st_read("~/Desktop/AGGREGATE-social/Scocial_ES_Mean_vector.shp") -> social_es_mean_shp

#try <- raster("~/Desktop/AGGREGATE-social/Social_ES_Sum_Rescale.tif")
#mapview(try)

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
#economic_es_mean_clip1_resam_raster <- raster(economic_es_mean_clip1_resam_filepath)
#mapview(social_es_mean_raster, query.digits = 1)
#maxValue(social_prioritized_raster)
ecological_prioritized_raster<- raster::raster(ecological_prioritized_filepath)
economic_prioritized_raster <- raster::raster(economic_prioritized_filepath)
social_prioritized_raster <- raster::raster(social_prioritized_filepath)
equal_weight_raster <- raster(equal_weight_filepath)
#mapview(economic_prioritized_raster)

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

minValue(economic_priori_hotspot_raster)
maxValue(social_prioritized_raster)
mapview(fran_adm1_shp)
mapview(economic_prioritized_raster)

# raster::crop(fran_adm1_shp, economic_prioritized_raster)
# y <- extent(fran_adm1_shp)
# satellite::crop(economic_prioritized_raster, y)
# crs(economic_prioritized_raster)
# st_crs(fran_adm1_shp)
# transform <- st_transform(fran_adm1_shp, crs = crs(economic_prioritized_raster))
# transform
# mapview(transform)
# economic_priori_crop <- crop(economic_prioritized_raster, transform)
# mapview(economic_priori_crop)

### normalize economic prioritized raster
norm_economic_prioritized_raster <- (economic_prioritized_raster - minValue(economic_prioritized_raster))/ (maxValue(economic_prioritized_raster) - minValue(economic_prioritized_raster)) 
mapview(economic_prioritized_raster)
#maxValue(norm_economic_prioritized_raster)
#mapview(norm_economic_prioritized_raster)
#minValue(economic_prioritized_raster)

## function to normalize rasters
norm_function <- function(x){
  (x - minValue(x)) / (maxValue(x)-minValue(x))
}

norm_function(economic_es_mean_clip1_raster) -> norm_economic_es_mean_clip1_raster
norm_function(social_es_mean_raster) -> norm_social_es_mean_raster
norm_function(ecological_es_mean_raster) -> norm_ecological_es_mean_raster
  
mapview(norm_ecological_es_mean_raster)  
mapview(ecological_prioritized_raster)
### normalize ecological prioritized raster
norm_ecological_prioritized_raster <- norm_function(ecological_prioritized_raster)

# normalize social prioritized raster
norm_social_prioritized_raster <- norm_function(social_prioritized_raster)

mapview(norm_social_es_mean_raster)
mapview(norm_social_prioritized_raster)
mapview(norm_ecological_prioritized_raster)
minValue(norm_ecological_prioritized_raster)
minValue(norm_economic_es_mean_clip1_raster)
maxValue(norm_social_es_mean_raster)

mapview(ecological_es_mean_raster)

#### calculate moran I for spatial autocorrelation
# positive values shows similar values are closely located 
# negative values shows dissimilar values are grouped together
# values ranges from 1 to -1 with 1 being perfect clustering and 0 being perfect randomness
economic_moran <- Moran(norm_economic_es_mean_clip1_raster)
ecological_moran <- Moran(ecological_es_mean_raster)
social_moran <- Moran(norm_social_es_mean_raster)

# MoranLocal(economic_es_mean_clip1_raster)
# Geary(economic_es_mean_clip1_raster)

freqDT(economic_es_mean_clip1_raster)
raster::corLocal(norm_economic_es_mean_clip1_raster, ecological_es_mean_raster)
compareRaster(norm_economic_es_mean_clip1_raster, norm_social_es_mean_raster)

raster::layerStats(economic_es_mean_clip1_raster, stat = 'pearson')
brick(norm_economic_es_mean_clip1_raster, ecological_es_mean_raster)
str(norm_economic_es_mean_clip1_raster)
str(ecological_es_mean_raster)

#### set rasters to have equal extent, crs, dimension to enable stacking and other analysis
## set attributes in reference to norm_ecological_es_mean_raster
ext_used <- extent(norm_ecological_es_mean_raster)
crs_used <- crs(norm_ecological_es_mean_raster)
res_used <- res(norm_ecological_es_mean_raster)
nrow_used <- dim(norm_ecological_es_mean_raster)[1]
ncol_used <- dim(norm_ecological_es_mean_raster)[2]

# set_norm_economic_es <- raster(norm_economic_es_mean_clip1_raster, values = TRUE, ext = ext_used,
#                                crs = crs_used, nrows = nrow_used, ncols = ncol_used)
# 
# try_economic <- raster(vals=values(norm_economic_es_mean_clip1_raster),ext=extent(norm_ecological_es_mean_raster),
#               crs=crs(norm_ecological_es_mean_raster),
#               nrows=dim(norm_ecological_es_mean_raster)[1],ncols=dim(norm_ecological_es_mean_raster)[2])

align_raster_function <- function(x, y){
  `extent<-`(x, extent(y)) -> x
  raster::nrow(x) <- raster::nrow(y)
  raster::ncol(x) <- raster::ncol(y)
}

##### align extent of norm_social_es_mean_raster with norm_ecological_es_mean_raster
norm_social_es_mean_raster <- alignExtent(extent = extent(norm_ecological_es_mean_raster), object = norm_social_es_mean_raster)

#norm_social_es_mean_raster <-  align_raster_function(norm_social_es_mean_raster, norm_ecological_es_mean_raster)

#align_raster_function(norm_economic_es_mean_clip1_raster)

`extent<-`(norm_economic_es_mean_clip1_raster, extent(norm_ecological_es_mean_raster)) -> norm_economic_es_mean_clip1_raster
#raster::`nrow<-`(norm_economic_es_mean_clip1_raster, nrow(norm_ecological_es_mean_raster)) -> norm_economic_es_mean_clip1_raster
raster::nrow(norm_economic_es_mean_clip1_raster) <- raster::nrow(norm_ecological_es_mean_raster)
raster::ncol(norm_economic_es_mean_clip1_raster) <- raster::ncol(norm_ecological_es_mean_raster)


`extent<-`(norm_social_es_mean_raster, extent(norm_ecological_es_mean_raster)) -> norm_social_es_mean_raster
raster::nrow(norm_social_es_mean_raster) <- raster::nrow(norm_ecological_es_mean_raster)
raster::ncol(norm_social_es_mean_raster) <- raster::ncol(norm_ecological_es_mean_raster)

######### set extent of norm_social_es_mean to norm_ecological_es_mean
norm_social_es_mean_raster <- setExtent(x = norm_social_es_mean_raster, ext = extent(norm_ecological_es_mean_raster))

norm_social_es_mean_raster_try <- raster(vals = values(norm_social_es_mean_raster), nrows = raster::nrow(norm_ecological_es_mean_raster),
                                     ncols = raster::ncol(norm_ecological_es_mean_raster), 
                                     resolution = res(norm_ecological_es_mean_raster))

res(norm_social_es_mean_raster)<-res_used
#res(norm_social_es_mean_raster, res_used) -> try_soc_norm_res
trysoc <- raster(social_es_mean_filepath, nrows = nrow_used, ext = ext_used)
`res<-`(trysoc, res_used)
raster::`nrow<-`(trysoc, nrow_used)
norm_social_es_mean_raster<- raster::setValues(trysoc)

# norm_social_es_mean_raster <- raster::`extent<-`(norm_social_es_mean_raster, trial_extent)
# 
# trial_extent<- extent(norm_ecological_es_mean_raster)

dim(norm_ecological_es_mean_raster)
dim(trysoc) <- dim(norm_ecological_es_mean_raster)

align_nrows <- function(x, y){
  raster::nrow(x) <- raster::nrow(y)
}

align_nrows(norm_social_es_mean_raster, norm_ecological_es_mean_raster)



all_es_mean_stack <- raster::stack(norm_economic_es_mean_clip1_raster, norm_ecological_es_mean_raster, norm_social_es_mean_raster)
raster::corLocal(all_es_mean_stack)
brick_economic_ecological<- raster::brick(norm_economic_es_mean_clip1_raster, norm_ecological_es_mean_raster)
#remove(norm_social_es_mean_raster)

alignExtent(ext_used, norm_economic_es_mean_clip1_raster)

tryalign <- alignExtent(extent(norm_ecological_es_mean_raster), norm_economic_es_mean_clip1_raster)
raster::corLocal(brick_economic_ecological)

raster::`extent<-`(norm_economic_es_mean_clip1_raster, ext_used)
norm_economic_es_mean_clip1_raster

terra::barplot(norm_economic_es_mean_clip1_raster)
terra::autocor(economic_es_mean_clip1_raster)
rastcombine <- c(norm_ecological_es_mean_raster, norm_economic_es_mean_clip1_raster)
typeof(rastcombine)
terra::compareGeom(norm_economic_es_mean_clip1_raster, norm_ecological_es_mean_raster)
rast(norm_economic_es_mean_clip1_raster)
rast(norm_ecological_es_mean_raster)

terra::density(norm_economic_es_mean_clip1_raster)
terra::density(norm_ecological_es_mean_raster)
terra::depth(rast(norm_ecological_es_mean_raster))
terra::ext(norm_ecological_es_mean_raster)
terra::pairs(rast(brick_economic_ecological), hist = T)
terra::persp(norm_ecological_es_mean_raster)


r22 <- raster(vals=values(r2),ext=extent(r1),crs=crs(r1),
              nrows=dim(r1)[1],ncols=dim(r1)[2])
