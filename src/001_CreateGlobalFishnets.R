library(raster)
library(tidyverse)
library(sf)
library(rgdal)
library(rgeos)
library(fasterize)
library(velox)
library(assertthat)
library(maptools)
require(gdalUtils)
source("src/000_ConvenienceFunctions.R")
raster::rasterOptions(progress = 'text')
raster::beginCluster(5)

# -------------------------------------------- #
#s.crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs" # World Mollweide
s.crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
latlong.crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" # Latitude longitude
res <- 1000                  # Source resolution, meters
tempdir <- '/media/martin/data/'
aggregate_source <- TRUE    # Should coarser grids be aggreated or created a new?
use_gadm <- TRUE            # Use GADM instead of natural earth
# -------------------------------------------- #
if(!dir.exists("data")) dir.create("data")

#### Create global fishnet temples ####
global_extent <- c(xmin = -179.999999, xmax = 180, ymin = -90, ymax = 90) # slightly smaller x because of dateline issues
# First create a global 1deg lat-long polygon
temp_latlong <- raster(xmn = global_extent[1] , xmx = global_extent[2],
                       ymn = global_extent[3], ymx = global_extent[4])
temp_latlong_fname <- paste0(tempdir,"latlong.tif")
writeGeoTiff(temp_latlong,temp_latlong_fname)

# For mollweide projection
temp_moll_fname <- paste0(tempdir,"moll.tif")
# Project to get the world extent to the Mollweide projection
gdalwarp(srcfile = temp_latlong_fname,
         dstfile = temp_moll_fname,
         s_srs = latlong.crs,
         t_srs = s.crs,
         r = "near",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         multi = TRUE,
         overwrite=TRUE,verbose=TRUE
)
global_extent_moll <- extent(raster(temp_moll_fname))

# Now to create the global template raster appropriate at full latitude - longitudes
global_raster <- raster(ext = global_extent_moll,
                        crs = s.crs,
                        resolution = res,
                        vals = 0
)
# Clean up
file.remove(temp_latlong_fname); file.remove(temp_moll_fname)
# --- #
if(use_gadm){
  world <- sf::read_sf('/media/martin/data/raw/gadm36_shp') %>% 
    dplyr::filter(NAME_0 != 'Antarctica') %>% 
    st_cast("MULTIPOLYGON") %>% 
    st_transform(crs = s.crs) %>% 
    mutate( dummy = 1 )
} else {
  library(rnaturalearth);library(rnaturalearthdata)
  # Continent factorization
  # - Only for natural earth data
  continents = c("Africa" = 1, "Antarctica" = 2, "Asia" = 3, "Europe" = 4, "North America" = 5,"Oceania" = 6, "Seven seas (open ocean)" = 7, "South America" = 8)
  # Get the NaturalEarthdata
  world <- ne_countries(scale = "large", returnclass = "sf") %>% 
    st_crop(c(xmin = -180, xmax=180,ymin= -90, ymax=90)) %>% 
    st_cast("MULTIPOLYGON") %>% 
    st_transform(crs = s.crs) %>% 
    mutate(continent_id = factor(continent,levels = names(continents),labels = as.vector(continents) )) %>% 
    mutate( dummy = 1 )
}
gc()

# ------------------- #
# Split the raster into 4 parts
dir.create( paste0(tempdir, 'temp_splits') )
writeGeoTiff(global_raster, paste0(tempdir, 'temp_splits/layer.tif') )
ss <- split_raster(paste0(tempdir, 'temp_splits/layer.tif'),
                   outpath = paste0(tempdir, 'temp_splits'),
                   2)

for(r in 1:nrow(ss)){
  print(r)
  # Use velox to keep small islands in
  #vx <- velox( raster(as.character( ss$file )[r]) )
  #vx$rasterize(spdf = world, field="dummy", background = NA, small = TRUE)
  #global_raster <- vx$as.RasterLayer();rm(vx)
  ref <- raster(as.character( ss$file )[r])
  global_raster <- fasterize::fasterize(sf = world,raster = ref)
  writeGeoTiff(global_raster, fname = paste0(tempdir, 'temp_splits/world_',r,'.tif')) # Save outputs and clear 
  rm(global_raster)
}
# Mosaic 
gdalUtils::gdalbuildvrt(gdalfile = list.files(paste0(tempdir, 'temp_splits'),"world",full.names = T),
                          output.vrt = paste0(tempdir, 'temp_splits/out.vrt'),
                          separate = FALSE)
global_raster <- raster(paste0(tempdir, 'temp_splits/out.vrt'))

writeGeoTiff(global_raster,fname = "data/globalgrid_mollweide_1km.tif") # Save outputs and clear 
# Load and reproject all layers to WGS84 and save again
# 1km
gdalwarp(srcfile = "data/globalgrid_mollweide_1km.tif",
         dstfile = "data/globalgrid_wgs84_1km.tif",
         s_srs = s.crs,
         t_srs = latlong.crs,
         r = "near",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         output_Raster=TRUE, multi = TRUE,
         overwrite=TRUE,verbose=TRUE)
lapply(list.files(paste0(tempdir, 'temp_splits')), file.remove) # clean up
unlink(paste0(tempdir, 'temp_splits'))

# ----------------------------- #
# Water coverage data from Copernicus 
# Used the permanent water occurence layer and subsetted to grid cells with at over 90% water coverage at 30m resolution
# This was then aggregated to 10km resolution by calculating the fractional value within each 10km grid cell
water <- raster::raster('/media/martin/data/Water_occurence_wgs84_1km.tif')
water <- raster::flip(water,'y')
# Correct caspian sea
casp <- sf::read_sf('/media/martin/data/CaspianSea.gpkg') %>% mutate(class = 1)
casp <- fasterize::fasterize(casp, water, field = "class")
water <- raster::mask(water, casp, maskvalue = 1, updatevalue = 0) # Mask out

# Now load and align global 1km grid
globalgrid <- raster('data/globalgrid_wgs84_1km.tif')
# Align water to global grid at 1km
water <- alignRasters(water, globalgrid,method = 'bilinear',func = mean, cl = T)
water <- raster::clamp(water,lower = 0, upper = 1000)
water <- raster::mask(water,globalgrid)
# Get the grid
#globalgrid <- globalgrid * 1000
globalgrid[globalgrid==1] <- 0
# Add them from each other to get the fractional cover at 1km WGS84
globalgrid <- (globalgrid + water)
writeGeoTiff(globalgrid, '/media/martin/data/Land_surface_wgs84_1km.tif')

# And for fractional cover at coarser resolution
land <- raster::raster('/media/martin/data/Land_surface_wgs84_1km.tif')
# Project to mollweide
projectMollWeide(land,'/media/martin/data/Land_surface_mollweide_1km.tif')
land <- raster('/media/martin/data/Land_surface_mollweide_1km.tif') # ~454.3084m resolution
# Align this to the 1km grid we ha
# Aggregate once and then align
land <- raster::aggregate(land, fact = 2, na.rm = TRUE, fun = mean)
# Average aggregate to 1km
gdalUtils::gdalwarp(srcfile = '/media/martin/data/Land_surface_mollweide_1km.tif',
                    dstfile = '/media/martin/data/globalgrid_mollweide_1km.tif',
                    s_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs',
                    t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs',
                    tr = c(1000,1000),
                    r = "average",
                    co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
                    multi = TRUE,
                    q = FALSE,overwrite=TRUE,verbose=TRUE
                    )
  
gc()
# ------------------------------------- #
# Now create aggregated layers at 10 and 50km resolution
beginCluster(7)
land <- raster('/media/martin/data/globalgrid_mollweide_1km.tif')
land <- setMinMax(land)

# Now aggregate
land[is.na(land)] <- 0
land5km_frac <- raster::aggregate(land, fact=5, fun=mean, expand=FALSE, na.rm=TRUE)
land5km_frac[land5km_frac==0] <- NA # Mask out fully water covered sites
writeGeoTiff(land5km_frac,fname = "/media/martin/data/globalgrid_mollweide_5km.tif") # Save outputs and clear 
# Then for 10km resolution
land10km_frac <- raster::aggregate(land, fact=10,expand=FALSE, na.rm=TRUE, fun = mean)
land10km_frac[land10km_frac==0] <- NA # Mask out fully water covered sites
writeGeoTiff(land10km_frac,fname = '/media/martin/data/globalgrid_mollweide_10km.tif') # Save outputs and clear 
# And for 50km
land50km_frac <- raster::aggregate(land, fact=50,expand=FALSE, na.rm=TRUE, fun = mean)
land50km_frac[land50km_frac==0] <- NA # Mask out fully water covered sites
writeGeoTiff(land50km_frac,fname = '/media/martin/data/globalgrid_mollweide_50km.tif') # Save outputs and clear 

# ----------------------- #
myLog('Now project all layers')
# 1km
gdalwarp(srcfile = "/media/martin/data/globalgrid_mollweide_1km.tif",
         dstfile = "/media/martin/data/globalgrid_wgs84_1km.tif",
         s_srs = s.crs,
         t_srs = latlong.crs,
         r = "bilinear",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         output_Raster=TRUE, multi = TRUE,
         overwrite=TRUE,verbose=TRUE)

# 5km
gdalwarp(srcfile = "/media/martin/data/globalgrid_mollweide_5km.tif",
         dstfile = "/media/martin/data/globalgrid_wgs84_5km.tif",
         s_srs = s.crs,
         t_srs = latlong.crs,
         r = "bilinear",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         output_Raster=TRUE, multi = TRUE,
         overwrite=TRUE,verbose=TRUE)

# 10km
gdalwarp(srcfile = "/media/martin/data/globalgrid_mollweide_10km.tif",
         dstfile = "/media/martin/data/globalgrid_wgs84_10km.tif",
         s_srs = s.crs,
         t_srs = latlong.crs,
         r = "bilinear",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         output_Raster=TRUE, multi = TRUE,
         overwrite=TRUE,verbose=TRUE)

# 50km
gdalwarp(srcfile = "/media/martin/data/globalgrid_mollweide_50km.tif",
         dstfile = "/media/martin/data/globalgrid_wgs84_50km.tif",
         s_srs = s.crs,
         t_srs = latlong.crs,
         r = "bilinear",
         co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
         output_Raster=TRUE, multi = TRUE,
         overwrite=TRUE,verbose=TRUE)

myLog('Created and projected global fishnet grids')
stop('DONE')
