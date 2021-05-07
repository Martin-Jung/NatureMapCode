# ------------------------------------------------------- #
# Prepare water data for subsequent analyses
library(raster)
library(fasterize)
library(gdalUtils)
library(rgdal)
library(velox)
library(sf)
library(stringr)
library(assertthat)
library(doParallel)
library(data.table)
library(tidyverse)
source("src/000_ConvenienceFunctions.R")
rasterOptions(progress = 'text')

data_path <- "data/"

# The path to the most up-to-date water data
water_path = "/media/martin/data/raw/Water"

target_resolution = "1km"; target_resolution_number <- 1
projection = "mollweide"

output_path = "/media/martin/data/features_esh/water"
if(!dir.exists(output_path)) dir.create(output_path)

# --------------------------------- # 
#### Realized water data ####
# Relative realised water provisioning services index
# Unit of 0 - 1
# Resolution is in ~10km

# Reproject
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","10km",".tif") )
fname_moll = paste0(output_path,"/Water_realised_","10km","_",projection,".tif") # Output name

gdalwarp(srcfile = paste0(water_path,"/realised_water_value.asc"),
         dstfile = fname_moll,
         r = "bilinear",
         s_srs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
         t_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",
         tr = raster::res(globalgrid),
         multi = TRUE,co=c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9")
)

# -------- Aggregate --------- #
water <- raster(fname_moll) # Security copy
water <- setMinMax(water)
water[is.na(water)] <- 0

water50km <- raster::aggregate(water, fact = 5,na.rm = TRUE, fun = mean)
fname_moll = paste0(output_path,"/Water_realised_","50km","_",projection,".tif") # Output name
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","50km",".tif") )
water50km <- alignRasters(water50km, globalgrid)
writeGeoTiff(water50km, fname_moll,"FLT4S") # Save outputs and clear 

assert_that(
  file.exists(fname_moll),
  cellStats(water50km,'max') == 1,
  compareRaster(water50km,globalgrid)
)

#### Potential clean water provisioning per basin ####
# Potential clean water provision normalized by basin
# Resolution is in ~10km

# Reproject
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","10km",".tif") )
fname_moll = paste0(output_path,"/cleanwaternorm_pot_","10km","_",projection,".tif") # Output name

gdalwarp(srcfile = paste0(water_path,"/Clean_water_norm_by_basin.tif"),
         dstfile = fname_moll,
         r = "bilinear",
         s_srs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
         t_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",
         tr = raster::res(globalgrid),
         multi = TRUE,co=c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE","ZLEVEL=9")
)
water <- raster(fname_moll)
water <- alignRasters(water, globalgrid)
writeGeoTiff(water, fname_moll, "FLT4S") # Save outputs and clear 

# -------- Aggregate --------- #
water <- raster(fname_moll) # Security copy
water <- setMinMax(water)
#water[water == 0] <- NA

water50km <- raster::aggregate(water, fact = 5,na.rm = TRUE, fun = mean)
fname_moll = paste0(output_path,"/cleanwaternorm_pot_","50km","_",projection,".tif") # Output name
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","50km",".tif") )
water50km <- alignRasters(water50km, globalgrid)
water50km <- raster::mask(water50km, globalgrid)
writeGeoTiff(water50km, fname_moll,"FLT4S") # Save outputs and clear 

assert_that(
  file.exists(fname_moll),
  cellStats(water50km,'max') <= 1,
  compareRaster(water50km,globalgrid)
)
