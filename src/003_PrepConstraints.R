# Purpose of this script is to prepare the constraints 
# for the global prioritization data

# Parameters
wdpa_data <- "/media/martin/data/raw/WDPA_Apr2019_Public.gpkg"
wdpa_point_data <- "/media/martin/data/raw/WDPA_Apr2019_Public_points.gpkg" # WDPA point dataset

# Ecoregion path
ecoregions_path <- "/media/martin/data/raw/Ecoregions2017/Ecoregions2017.shp"

# Artifical mask
art <- '/mnt/hdrive/IUCNHabitatmapping/ForExport/lvl1_frac_1km_ver003/iucn_habitatclassification_fraction_lvl1__1400_Artificial - Terrestrial__ver003.tif' # at 1km resolution

# Target resolution identifier
target_resolution <- "50km"
target_resolution_number <- 50
projection = "mollweide"
wdpa_includepoints = FALSE # Should the buffered point data be included in the calculation
wdpa_olderversion = c(NA,'2010','2014')[1] # Subset the WDPA data to STATUS_YEAR before/equal 2014

output_path <- "/media/martin/data/constraints"

# Output paths
if(!dir.exists(output_path)) dir.create(output_path)

# ------------------------- #
# Start processing #
# ------------------------- #
library(rgdal)
library(sf)
library(raster)
library(velox)
library(lwgeom)
library(gdalUtils)
library(fasterize)
library(assertthat)
source('src/000_ConvenienceFunctions.R')
# Security checks
stopifnot(
  assert_that(
    file.exists(wdpa_data),
    file.exists(wdpa_point_data),
    file.exists(kba_data),
    file.exists(ecoregions_path)
  )
)

myLog('Load target resolution grid and input files')

# The global grid
global_raster1km <- raster(paste0('data/globalgrid_',projection,'_1km.tif'))
global_raster1km[global_raster1km>1] <- 1
global_raster <- raster(paste0('data/globalgrid_',projection,'_',target_resolution,'.tif'))
global_raster[global_raster > 1] <- 1

# Matchup table for dummy variable
dummy_table <- tibble::tibble(
  IUCN_CAT = c("Ia","Ib","II","III","IV","V","VI","Not Applicable","Not Reported","Not Assigned"," "),
  dummy = c(1,2,3,4,5,6,7,11,12,13,12)
)

# ---------------------------------------------- #
# Should we do a status year subset?
if(!is.na(wdpa_olderversion)){
  myLog('Loading older WDPA version...')
  # Load older versions of WDPA
  wp <- paste0("/media/martin/data/raw/WDPA_",wdpa_olderversion,"_poly.gpkg")
  stopifnot( assert_that(file.exists(wp)) )
    
  wdpa <-  sf::st_read(wp) %>% 
    dplyr::filter(STATUS %notin% c('Proposed', 'Not Reported') ) %>% # Remove proposed and not reported PA's
    dplyr::filter(DESIG_ENG != 'UNESCO-MAB Biosphere Reserve') %>% # Remove MAB sites
    dplyr::filter(ISO3 %notin% c('CHN','EST')) %>% # Remove china and Estonia
    st_cast("MULTIPOLYGON") 

} else {
  myLog('Loading current WDPA version...')
  # Load public WDPA polygon data
  wdpa <- sf::st_read(wdpa_data) %>% 
    dplyr::filter(STATUS %notin% c('Proposed', 'Not Reported') ) %>% # Remove proposed and not reported PA's
    dplyr::filter(DESIG_ENG != 'UNESCO-MAB Biosphere Reserve') %>% # Remove MAB sites
    st_cast("MULTIPOLYGON") 
} 

# Create a new dummy variable with the categories in them
wdpa <- dplyr::left_join(wdpa, dummy_table, by = "IUCN_CAT") %>% 
  dplyr::select(STATUS_YR,dummy)

if(wdpa_includepoints){
  wdpa_points <- sf::st_read(wdpa_point_data)
  # Reproject then buffer using meter
  x = wdpa_points %>% st_transform(crs = proj4string(global_raster)) %>% 
    dplyr::filter(DESIG_ENG != 'UNESCO-MAB Biosphere Reserve') %>% # Remove MAB sites
    # Remove points without reported area
    dplyr::filter(REP_AREA > 0) %>% 
    # Construct the buffered area field to square meter assuming the area is a 'circle'
    dplyr::mutate(buffer = sqrt((REP_AREA * 1e6)/pi) )
  # Now buffer and transform back
  wdpa_points <- sf::st_buffer(x,dist = x$buffer) %>% 
    st_cast("MULTIPOLYGON") %>% 
    st_transform(crs = st_crs(wdpa)$proj4string)
  
  # Make dummy columns
  wdpa_points <- wdpa_points %>% dplyr::mutate(dummy = dummy_table$dummy[which(dummy_table$IUCN_CAT=="Not Reported")]
                                               )
  # Select only columns needed
  wdpa_points <- wdpa_points %>% dplyr::select(STATUS_YR,dummy)
  
  # Append
  wdpa <- sf:::rbind.sf(wdpa,wdpa_points)
  rm(wdpa_points,x)
}

# Make valid
wdpa <- lwgeom::st_make_valid(wdpa)
wdpa <- wdpa %>% st_cast("MULTIPOLYGON") 

#box = c(xmin = -179.9999, ymin = -90, xmax = 180, ymax = 90)
#st_crop(box) %>% st_cast("MULTIPOLYGON")

if(projection!= "wgs84"){
  wdpa <- st_transform(wdpa,crs = proj4string(global_raster) ) #options = "WRAPDATELINE=NO")
}
# Check plot  
plot(wdpa['dummy'])

stopifnot(
  assert_that( length( which(is.na(wdpa$dummy)) ) == 0 )
)
# ------------------------------ #
gc(verbose = TRUE);myLog('Now rasterizing to grid')

# First rasterize to finest resolution grid (e.g. 1km)
#vx <- velox(global_raster1km)
#vx2 <- vx$copy()
#vx2$rasterize(spdf=wdpa , field="dummy", background=NA, small = TRUE)
#wdpa_raster <- vx2$as.RasterLayer();rm(vx2);gc()
# gdalUtils::gdal_rasterize(src_datasource = "/media/martin/data/wdpa_temp.shp",
#                           dst_filename = paste0(output_path,'/wdpa_april2019_','1km',
#                                                 ifelse(projection=="wgs84","",paste0("_",projection)),
#                                                 '.tif'),
#                           te = as.ve  ctor(extent(global_raster1km)),
#                           ts = res(global_raster1km),
#                           ot = "INT1U",
#                           a = "dummy",
#                           co = c("COMPRESS=DEFLATE","ZLEVEL=9")
#                           )
wdpa_raster <- fasterize::fasterize(wdpa, global_raster1km,field ="dummy", fun = "first")
writeGeoTiff(wdpa_raster, 
             paste0(output_path,'/wdpa_',
                                 ifelse(is.na(wdpa_olderversion),'april2019',wdpa_olderversion),'_','1km',
                                 ifelse(projection=="wgs84","",paste0("_",projection)),
                                 '.tif')
             ) # Save outputs and clear 

# Then aggregate (1) modal  
wdpa_raster_copy <- wdpa_raster # Security copy
wdpa_raster[is.na(wdpa_raster)] <- 0

# And aggregate (2) as fraction of cell
# sum all values that equal 1, and divide by the total number of pixels that are being aggregated (10*10)
wdpa_raster[wdpa_raster>=1] <- 1
wdpa_raster10km <- raster::aggregate(wdpa_raster, 
                                     fact = target_resolution_number, 
                                     na.rm = TRUE,
                                     fun = function(x, ...)
                                     { (sum(x == 1, ...)/ (target_resolution_number*target_resolution_number))}
)
wdpa_raster10km[wdpa_raster10km == 0] <- NA # Replace all non-PA with NA again
writeGeoTiff(wdpa_raster10km, 
             paste0(output_path,'/wdpa_',
                    ifelse(is.na(wdpa_olderversion),'april2019',wdpa_olderversion),'_',target_resolution,
                    ifelse(projection=="wgs84","",paste0("_",projection)),
                                     '_fraction.tif'),
             'FLT4S') # Save outputs and clear 
rm(wdpa_raster,wdpa_raster10km);gc() # Cleanup


#### Biomes and ecoregions ####
# Create rasterized version of global biomes and ecoregions
# Load ecoregions
ecoregions <- sf::st_read(ecoregions_path)

if(projection!= "wgs84"){
  ecoregions <- st_transform(ecoregions,crs = proj4string(global_raster))
}

# --- Biome --- #
# First rasterize to finest resolution grid (e.g. 1km)
#vx <- velox(global_raster1km)
#vx2 <- vx$copy()
#vx2$rasterize(spdf=ecoregions , field="BIOME_NUM", background=NA, small = TRUE)
#biome_raster <- vx2$as.RasterLayer();rm(vx2);gc()
biome_raster <- fasterize(ecoregions, global_raster1km,field ="BIOME_NUM",fun = "first")
writeGeoTiff(biome_raster, paste0(output_path,'/biomes_','1km',
                                  ifelse(projection=="wgs84","",paste0("_",projection)),
                                  '.tif')
             ) # Save outputs and clear 

# Then aggregate (1) modal
biome_raster_copy <- biome_raster # Security copy
biome_raster10km <- raster::aggregate(biome_raster, fact = target_resolution_number,
                                      na.rm = TRUE, fun = raster::modal)
# Crop and align to globalgrid
biome_raster10km <- raster::crop(biome_raster10km,global_raster)
biome_raster10km <- alignRasters(biome_raster10km,template = global_raster,method = 'ngb',func = raster::modal,cl = T)
writeGeoTiff(biome_raster10km, paste0(output_path,'/biomes_',target_resolution,
                                      ifelse(projection=="wgs84","",paste0("_",projection)),
                                      '_modal.tif')
             ) # Save outputs and clear 
rm(biome_raster10km);gc()

# --- Ecoregion --- #
#vx <- velox(global_raster1km)
#vx2 <- vx$copy()
#vx2$rasterize(spdf=ecoregions , field="ECO_ID", background=NA, small = TRUE)
#eco_raster <- vx2$as.RasterLayer();rm(vx2);gc()
eco_raster <- fasterize(ecoregions, global_raster1km,field ="ECO_ID", fun = "first")
writeGeoTiff(eco_raster, paste0(output_path,'/ecoregions_','1km',
                                ifelse(projection=="wgs84","",paste0("_",projection)),
                                '.tif')
             ) # Save outputs and clear 

# Then aggregate (1) modal
eco_raster_copy <- eco_raster # Security copy
eco_raster10km <- raster::aggregate(eco_raster, fact = target_resolution_number,
                                    na.rm = TRUE, fun = raster::modal)
eco_raster10km <- raster::crop(eco_raster10km,global_raster)
eco_raster10km <- alignRasters(eco_raster10km,template = global_raster, method = 'ngb', func = raster::modal, cl = T)
writeGeoTiff(eco_raster10km, paste0(output_path,'/ecoregions_',target_resolution,
                                    ifelse(projection=="wgs84","",paste0("_",projection)),
                                    '_modal.tif')
             ) # Save outputs and clear 
rm(eco_raster10km);gc()

## --- ##
#### POWO TDWG regions ####
pu_id10km <- raster('data/puid_10km.tif')
pu_id50km <- raster('data/puid_50km.tif')

library(fst)
powo_shape <- sf::read_sf("/mnt/hdrive/Otherstuff/PlantsMapping/wgsrpd-master/level2/level2.shp")  %>% 
  dplyr::select(LEVEL2_COD)
# Make valid
powo_shape <- lwgeom::st_make_valid(powo_shape) %>% st_cast("MULTIPOLYGON") 
# Crop to specific bbox because of dateline issues
box = c(xmin = -179.9999, ymin = -90, xmax = 180, ymax = 90)
powo_shape <- st_crop(powo_shape, box) %>% st_cast("MULTIPOLYGON")
powo_shape <- powo_shape %>% st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")


powo_shape_2 <- fasterize(powo_shape,pu_id10km,field = 'LEVEL2_COD')

powo_shape <- sf::read_sf("/mnt/hdrive/Otherstuff/PlantsMapping/wgsrpd-master/level3/level3.shp")  %>% 
  dplyr::select(LEVEL3_COD)
# Make valid
powo_shape <- lwgeom::st_make_valid(powo_shape) %>% st_cast("MULTIPOLYGON") 
# Crop to specific bbox because of dateline issues
box = c(xmin = -179.9999, ymin = -90, xmax = 180, ymax = 90)
powo_shape <- st_crop(powo_shape, box) %>% st_cast("MULTIPOLYGON")
powo_shape <- powo_shape %>% st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

powo_shape$LEVEL3_NUM <- as.numeric( as.factor(powo_shape$LEVEL3_COD) )
powo_shape_3 <- fasterize(powo_shape,pu_id10km,field = 'LEVEL3_NUM')

df <- data.frame(wgsrpd_lvl2 = values(powo_shape_2),wgsrpd_lvl3 = values(powo_shape_3), pu_id_10km = values(pu_id10km) )
df$wgsrpd_lvl2[which(df$pu_id_10km==0)] <- NA
df$wgsrpd_lvl3[which(df$pu_id_10km==0)] <- NA
df <- subset(df, complete.cases(df))
# Join in ids again 
df <- dplyr::left_join(df, powo_shape %>%  sf::st_drop_geometry() %>% as.data.frame(), by = c('wgsrpd_lvl2' = 'LEVEL3_NUM') )
df$wgsrpd_lvl3 <- df$LEVEL3_COD
df$LEVEL3_COD <- NULL

# Write out put
write.fst(df,'data/pu_wgsrpd_50km.fst')

powo_shape$pu_id_10km <- raster::extract(pu_id10km, powo_shape)
powo_shape$pu_id_50km <- raster::extract(pu_id50km, powo_shape)

# Write out put
write.fst(as.data.frame(powo_shape %>% sf::st_drop_geometry()  ),
          'data/pu_wgsrpd_10km.fst')


#### Artificial masking ####
# Create updated mask.
# Ensure that amount of land and total artifical land harmonizes

# Land surface area globally
land <- raster('/media/martin/data/Land_surface_wgs84_1km.tif')
# Artifical area globally
ras_art <- setMinMax( raster(art) ) # set from above
#ras_art[is.na(ras_art)] <- 0 # overwrite with zero
# Align them both
ras_art <- alignRasters(ras_art, land, method = 'bilinear', func = mean)
# Now get the minimum from both
#new <- raster::calc( stack(ras_art,land),fun = min, na.rm = TRUE )
new <- ras_art
# Finally mask with land area
land.mask <- land
land.mask[land.mask >0] <- 1
new <- raster::mask(new, land.mask)

assert_that(
  exists('new'),
  cellStats(new,'max') <= 1000
)

# For 10km
new_raster10km <- raster::aggregate(new, fact = 10,
                                         na.rm = TRUE, fun = mean)
writeGeoTiff(new_raster10km, paste0('/media/martin/data/Artifical_lvl1_10km_mollweide_11082020.tif'),"INT4S") # Save outputs and clear 
# And for 50km
new_raster50km <- raster::aggregate(new, fact = 50,
                                         na.rm = TRUE, fun = mean)
writeGeoTiff(new_raster50km, paste0('/media/martin/data/Artifical_lvl1_50km_mollweide_11082020.tif'),"INT4S") # Save outputs and clear 
