# ------------------------------------------------------- #
# Prepare carbon data for subsequent analyses
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

# The path to the most up-to-date carbon data
carbon_path = "/media/martin/data/raw/Carbon/"
carbon_path_geobene = "/mnt/pdrive/geobene2/Spatial_data/global/NatureMap/Biomass"

target_resolution = "10km"; target_resolution_number <- 10
projection = "mollweide"

#output_path = "/media/martin/data/features/carbon_agbc"
output_path = "/media/martin/data/features_esh/carbon_agbc"
if(!dir.exists(output_path)) dir.create(output_path)

# Load template raster
# --------------------------------- # 
#### SOC Data ####
# Carbon (merge beforehand in QGIS)
# Align with 1km global grid
globalgrid <- raster( paste0(data_path,"globalgrid_","wgs84","_","10km",".tif") )
ras <- raster("/media/martin/data/raw/Carbon/SOC_10km_v11.tif")
ras <- raster::flip(ras,direction = 'y')
ras <- alignRasters(ras, globalgrid, method = 'bilinear',func = raster::mean,cl = TRUE)

# Mask out using the terrestrial layer
m = globalgrid
m[m > 0] <- 1
ras <- raster::mask(ras, m)
ras <- raster::clamp(ras, lower = 0, upper = Inf)
writeGeoTiff(ras, "/media/martin/data/raw/Carbon/SOC_10km.tif", dt = 'FLT4S')

# Reproject
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","10km",".tif") )
fname <- paste0(carbon_path,"SOC_10km.tif")
fname_moll = paste0(output_path,"/Carbon_SOC_","10km","_",projection,".tif") # Output name
gdalwarp(srcfile = fname,dstfile = fname_moll,
         r = "bilinear",
         s_srs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
         t_srs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs",
         tr = raster::res(globalgrid),
         dstnodata = -9999,
         multi = TRUE,co=c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9")
)
# Align grid
carb <- raster(fname_moll)
carb <- alignRasters(carb, globalgrid, method = "bilinear",func = mean, cl = FALSE)
writeGeoTiff(carb, fname_moll,"FLT4S") # Save outputs and clear 
rm(globalgrid,carb)

# -------- Aggregate to 50km --------- #
carb <- raster(fname_moll) # Security copy
carb <- setMinMax(carb)
carb[is.na(carb)] <- 0
#carb10km <- raster::aggregate(carb, fact = 10,na.rm = TRUE, fun = mean)
#fname_moll = paste0(output_path,"/Carbon_SOC_","10km","_",projection,".tif") # Output name
#writeGeoTiff(carb10km, fname_moll,"FLT4S") # Save outputs and clear 

carb50km <- raster::aggregate(carb, fact = 5,na.rm = TRUE, fun = mean)
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","50km",".tif") )
carb50km <- alignRasters(carb50km, globalgrid, method = "bilinear",func = mean, cl = FALSE)
carb50km[carb50km==0] <- NA
fname_moll = paste0(output_path,"/Carbon_SOC_","50km","_",projection,".tif") # Output name
writeGeoTiff(carb50km, fname_moll, "FLT4S") # Save outputs and clear 

#### AGBC Data ####
# Below ground carbon
# Making a global ABGC layer in T/ha
ras_abgc <- raster(paste0(carbon_path_geobene,"/","AGBC_fmc.tif")) # Above and below ground biomass carbon
globalgrid <- raster( paste0(data_path,"globalgrid_","wgs84","_","1km",".tif") )

# Output name
fname = paste0(output_path,"/Carbon_AGBC_","1km","_","wgs84",".tif") # Output name
assert_that(!file.exists(fname))
gdalwarp(srcfile = paste0(carbon_path_geobene,"/","AGBC_fmc_1.tif"),dstfile = fname,
         r = "bilinear",
         tr = raster::res(globalgrid),
         multi = TRUE,
         co=c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9")
)
# Now load the 1km version and aggregate to 10km
ras_abgc1km <- raster(fname)
#ras_abgc1km[ras_abgc1km==0] <- NA
ras_abgc10km <- raster::aggregate(ras_abgc1km, fact = 10,na.rm = TRUE, fun = mean)
ras_abgc50km <- raster::aggregate(ras_abgc1km, fact = 50,na.rm = TRUE, fun = mean)
fname = paste0(output_path,"/Carbon_AGBC_","10km","_","wgs84",".tif") # Output name
writeGeoTiff(ras_abgc10km, fname, dt = 'FLT4S')
fname = paste0(output_path,"/Carbon_AGBC_","50km","_","wgs84",".tif") # Output name
writeGeoTiff(ras_abgc50km, fname, dt = 'FLT4S')

# Reproject to Mollweide and repeat aggregation
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","1km",".tif") )
(fname_moll = paste0(output_path,"/Carbon_AGBC_","1km","_",projection,".tif")) # Output name
ras <- raster(paste0(output_path,"/Carbon_AGBC_","1km","_","wgs84",".tif"))
carb <- projectRaster(ras, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
carb <- alignRasters(carb, globalgrid, method = "bilinear",func = mean)
carb <- raster::mask(carb,globalgrid)
writeGeoTiff(carb, fname_moll, dt = 'FLT4S')

# -------- Aggregate --------- #
#carb[carb==0] <- NA # Unnecessary, already 0
carb10km <- raster::aggregate(carb, fact = 10,na.rm = TRUE, fun = mean)
# Align with globalgrid 
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","10km",".tif") )
carb10km <- alignRasters(carb10km, globalgrid, method = 'bilinear',func = mean,cl = T)

fname_moll = paste0(output_path,"/Carbon_AGBC_","10km","_",projection,".tif") # Output name
carb10km[carb10km==0] <- NA
writeGeoTiff(carb10km, fname_moll,"FLT4S") # Save outputs and clear 

carb50km <- raster::aggregate(carb, fact = 50,na.rm = TRUE, fun = mean)
# Align with globalgrid 
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","50km",".tif") )
carb50km <- alignRasters(carb50km,globalgrid,method = 'bilinear',func = mean,cl = T)
carb50km[carb50km==0] <- NA
fname_moll = paste0(output_path,"/Carbon_AGBC_","50km","_",projection,".tif") # Output name
writeGeoTiff(carb50km, fname_moll,"FLT4S") # Save outputs and  clear 

# ---------------------------- #
#### BGBC Data ####
# Below ground carbon
# Making a global ABGC layer in T/ha
ras_abgc <- raster(paste0(carbon_path_geobene,"/","BGBC.tif")) # Above and below ground biomass carbon
globalgrid <- raster( paste0(data_path,"globalgrid_","wgs84","_","1km",".tif") )

# Output name
fname = paste0(output_path,"/Carbon_BGBC_","1km","_","wgs84",".tif") # Output name
assert_that(!file.exists(fname))
gdalwarp(srcfile = paste0(carbon_path_geobene,"/","BGBC.tif"),dstfile = fname,
         r = "bilinear",
         tr = raster::res(globalgrid),
         multi = TRUE,
         co=c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9")
)

# Reproject to Mollweide and repeat aggregation
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","1km",".tif") )
(fname_moll = paste0(output_path,"/Carbon_BGBC_","1km","_",projection,".tif")) # Output name
ras <- raster(paste0(output_path,"/Carbon_BGBC_","1km","_","wgs84",".tif"))
carb <- projectRaster(ras, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
carb <- alignRasters(carb, globalgrid, method = "bilinear",func = mean)
carb <- raster::mask(carb,globalgrid)
writeGeoTiff(carb, fname_moll, dt = 'FLT4S')

# -------- Aggregate --------- #
#carb[carb==0] <- NA # Unnecessary, already 0
carb10km <- raster::aggregate(carb, fact = 10,na.rm = TRUE, fun = mean)
# Align with globalgrid 
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","10km",".tif") )
carb10km <- alignRasters(carb10km, globalgrid, method = 'bilinear',func = mean,cl = T)

fname_moll = paste0(output_path,"/Carbon_BGBC_","10km","_",projection,".tif") # Output name
carb10km[carb10km==0] <- NA
writeGeoTiff(carb10km, fname_moll,"FLT4S") # Save outputs and clear 

carb50km <- raster::aggregate(carb, fact = 50,na.rm = TRUE, fun = mean)
# Align with globalgrid 
globalgrid <- raster( paste0(data_path,"globalgrid_",projection,"_","50km",".tif") )
carb50km <- alignRasters(carb50km,globalgrid,method = 'bilinear',func = mean,cl = T)
carb50km[carb50km==0] <- NA
fname_moll = paste0(output_path,"/Carbon_BGBC_","50km","_",projection,".tif") # Output name
writeGeoTiff(carb50km, fname_moll,"FLT4S") # Save outputs and clear 
