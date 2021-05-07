# ------------------------------------------------------- #
# The preparation script for all the models we refitted from BIEN + GBIF + INAT data  
library(raster)
library(fasterize)
library(rgdal)
library(gdalUtils)
library(NatureMapRangeModelling)
library(velox)
library(sf)
library(stringr)
library(assertthat)
library(doParallel)
library(data.table)
library(tidyverse)
source("src/000_ConvenienceFunctions.R")
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) # Chunking script
matchName2POWO <- function(plants, powo_checklist, powo_checklist_acc){
  assertthat::assert_that(
    is.data.frame(plants),
    is.data.frame(powo_checklist),
    is.data.frame(powo_checklist_acc),
    has_name(plants,'binomial'),
    has_name(powo_checklist_acc,"acc_full_name_without_family"),
    has_name(powo_checklist, "full_name_without_family")
  )
  # Convert to lower for matching
  plants$binomial <- str_to_lower(plants$binomial)
  # Select only POWO species in the list
  # First by getting the id per species
  plants <- dplyr::left_join(plants, powo_checklist_acc %>% dplyr::select(acc_full_name_without_family,accepted_db_id) %>% distinct() %>% 
                               dplyr::mutate(acc_full_name_without_family = str_to_lower(str_replace_all(acc_full_name_without_family," ","_"))),
                             by = c("binomial" = "acc_full_name_without_family") ) %>% 
    dplyr::rename('db_id' = 'accepted_db_id') 
  
  # Now also match the synonms
  # Now check with all the synonym's
  plants_syn <- dplyr::filter(plants,is.na(db_id)) %>% dplyr::select(-db_id) %>% 
    dplyr::left_join(., powo_checklist %>% dplyr::select(full_name_without_family,acc_full_name_without_family) %>% 
                       dplyr::mutate(acc_full_name_without_family = str_to_lower( str_replace_all(acc_full_name_without_family," ","_"))),
                     by = c("binomial" = "full_name_without_family")) %>% 
    drop_na(acc_full_name_without_family)  # Drop all remaining without accepted species name
  # Join in db_id of accepted species again
  plants_syn <- dplyr::left_join( plants_syn, powo_checklist_acc %>% dplyr::select(acc_full_name_without_family,accepted_db_id) %>% distinct() %>% 
                                    dplyr::rename('db_id' = 'accepted_db_id') %>% 
                                    dplyr::mutate(acc_full_name_without_family = str_to_lower( str_replace_all(acc_full_name_without_family," ","_"))),
                                  by = "acc_full_name_without_family") %>% 
    # Remove again all those that could not be matched
    drop_na(db_id)
  
  # Now append them both
  plants <- bind_rows(
    plants %>% dplyr::filter(!is.na(db_id)),
    plants_syn
  ) %>% dplyr::select(-acc_full_name_without_family)
  rm(plants_syn)
  if(nrow(plants)==0){ warning('No names could be matched!') }
  assert_that(!anyNA(plants$db_id))
  return(plants)
}
cores <- 7 # Number of cores and cluster_type for parallel processing 

# Global grid path
gridpath = "data/"
projection = "mollweide"

# POWO data
process_only_powo_species <- TRUE # Only process species that can be assigned to POWO taxonomy
powo_checklist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/checklist_species.txt') %>% 
  dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
powo_checklist_acc <- powo_checklist %>% dplyr::filter(taxon_status_description == 'Accepted')
powo_dist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/dist_NatureMap.txt') %>% 
  dplyr::filter(db_id %in% powo_checklist_acc$db_id )

# Path to processed files
newbien_path <- paste0("/media/martin/data/plant_sdms/SDMruns")

# Some error checks
stopifnot(
  assert_that(
    file.exists(newbien_path),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif"))
  )
)

# Fire up cores
if(file.exists('newBIEN_log.txt')) {file.remove('newBIEN_log.txt')}
cl <- makeCluster(cores, outfile = 'newBIEN_log.txt')
registerDoParallel(cores = cores)

# ------------------------------------------------------- #
# Make a global land area map out of the grid
globalgrid <- raster(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif"))
globalgrid <- setMinMax(globalgrid)
globalgrid[globalgrid>=1] <- 1
globalgrid50km <- raster(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif"))
# And the fraction cover
globalgrid_fraction <- raster(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")) / 1000
globalgrid_fraction50km <- raster(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif")) / 1000

# ------------------------------------------------------- #
# Load all the new layers
ll_ppm1 <- list.files(paste0(newbien_path,"/palearctic_biome/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE)
ll_ppm2 <- list.files(paste0(newbien_path,"/palearctic_inat/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE)
ll_ppm <- c( ll_ppm1[ which( basename(ll_ppm1) %notin% basename(ll_ppm2) ) ], ll_ppm2 ) # Any INAT added on top
ll_ppm <- c( ll_ppm, list.files(paste0(newbien_path,"/NZ_landcare/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_ppm <- c( ll_ppm, list.files(paste0(newbien_path,"/ian_southeastasia/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_ppm <- c( ll_ppm, list.files(paste0(newbien_path,"/african_endemics_inat/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_ppm <- c( ll_ppm, list.files(paste0(newbien_path,"/IM1/PPM/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
rm(ll_ppm1, ll_ppm2)
# RB
ll_rangebag1 <- list.files(paste0(newbien_path,"/palearctic_biome/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE)
ll_rangebag2 <- list.files(paste0(newbien_path,"/palearctic_inat/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE)
ll_rangebag <- c( ll_rangebag1[ which( basename(ll_rangebag1) %notin% basename(ll_rangebag2) ) ], ll_rangebag2 ) # Any INAT added on top
ll_rangebag <- c( ll_rangebag, list.files(paste0(newbien_path,"/NZ_landcare/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_rangebag <- c( ll_rangebag, list.files(paste0(newbien_path,"/ian_southeastasia/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_rangebag <- c( ll_rangebag, list.files(paste0(newbien_path,"/african_endemics_inat/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_rangebag <- c( ll_rangebag, list.files(paste0(newbien_path,"/IM1/RangeBag/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
rm(ll_rangebag1, ll_rangebag2)
ll_rangebag <- ll_rangebag[str_detect(basename(ll_rangebag),"X0.165")] # Grep only those that with 16.5% thresholding
# Points
ll_points1 <- list.files(paste0(newbien_path,"/palearctic_biome/Points/BinaryMaps"),"*.tif",recursive = T,full.names = TRUE)
ll_points2 <- list.files(paste0(newbien_path,"/palearctic_inat/Points/BinaryMaps"),"*.tif",recursive = T,full.names = TRUE)
ll_points <- c( ll_points1[ which( basename(ll_points1) %notin% basename(ll_points2) ) ], ll_points2 ) # Any INAT added on top
ll_points <- c( ll_points, list.files(paste0(newbien_path,"/NZ_landcare/Points/BinaryMaps"),"*.tif",recursive = T,full.names = TRUE) )
ll_points <- c( ll_points, list.files(paste0(newbien_path,"/ian_southeastasia/Points/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_points <- c( ll_points, list.files(paste0(newbien_path,"/african_endemics_inat/Points/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
ll_points <- c( ll_points, list.files(paste0(newbien_path,"/IM1/Points/BinaryMaps/"),"*.tif",recursive = T,full.names = TRUE) )
rm(ll_points1,ll_points2)

assert_that(
    length(ll_ppm) > 0,
    length(ll_rangebag) > 0,
    length(ll_points) > 0
  )

#### New PPM models ####
# Load all the PPM estimated ranges
ll <- ll_ppm
# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = str_split(basename(ll),"__",simplify = T)[,2], fname = as.character(ll) ,stringsAsFactors = FALSE) 
plants$binomial <- str_replace(plants$binomial,"\\.","_")

# Detect subspecies and correct subspecies to species
plants$subspecies <- Reduce("|", list(str_detect(plants$binomial,'subsp'), str_detect(plants$binomial,'var.')) )
plants_subspecies <- plants %>% dplyr::filter(subspecies == TRUE)
plants_subspecies$binomial <- paste(str_split(plants_subspecies$binomial,"_",simplify = T)[,1],str_split(plants_subspecies$binomial,"_",simplify = T)[,2])
# Join back
plants <- bind_rows( 
  plants %>% dplyr::filter(subspecies == FALSE),
  plants_subspecies
)
rm(plants_subspecies)

# And the number of duplicated species
species_duplicates <- names(which(table(plants$binomial)>1)) 
# Number of duplicates
myLog("Number of globally duplicated species: ",length( species_duplicates ))

if(process_only_powo_species){
  plants <- matchName2POWO(plants, powo_checklist, powo_checklist_acc)
  # Also limit powo dist to the target ids
  powo_dist <- powo_dist %>% dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}
data("tdwg")

# Process first the non-duplicated species!
plants_unique <- plants

# Output path
output_path_10km = paste0("/media/martin/data/features/newPlantPPM_","10km")
output_path_50km = paste0("/media/martin/data/features/newPlantPPM_","50km")
if(!dir.exists(output_path_10km)){ dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)){ dir.create(output_path_50km) }

# ------------------------ #
o <- foreach::foreach( row = 1:nrow(plants_unique),
                       .inorder = F,
                       .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                       .export = c("ll","globalgrid","plants_unique",
                                   "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                                   "output_path_10km","output_path_50km","powo_dist","tdwg",
                                   "alignRasters","writeGeoTiff")
) %dopar% {
  out_10km = paste0(output_path_10km,"/",plants_unique$binomial[row],".tif")
  out_50km = paste0(output_path_50km,"/",plants_unique$binomial[row],".tif")
  if((!file.exists(out_10km)) & (!file.exists(out_50km)) ){
    species_range <- raster(plants_unique$fname[row])
    if(is.infinite(cellStats(species_range,'max'))) { return(NULL)} # 33005 is empty?
    
    print(paste0("Preparing species ",plants_unique$binomial[row]))
    
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == plants_unique$db_id[row])
    
    # Reproject
    species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
    
    # Extend to resolution of global grid
    species_range <- raster::extend(species_range, globalgrid)
    
    # Align resolution
    species_range[species_range==0] <- NA
    species_range2 <- alignRasters(species_range,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
    
    if( raster::cellStats(species_range2,'max') != 1 ){
      # Bilinear cell interpolation seems to work in cases ngb fails
      species_range <- alignRasters(species_range,globalgrid, method = "bilinear",func = raster::modal, cl = FALSE)
      species_range[species_range2>0] <- 1
      rm(species_range2)
    } else { species_range <- species_range2; rm(species_range2) }
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid),
      raster::cellStats(species_range,'max') == 1
    )
    
    # ----- #
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    if(raster::cellStats(species_range,'max') == -Inf) { return(NULL)} # Filtering removed all habitat patches (unlikely)
    # ----- #
    
    # Also aggregate to 50km
    species_range50km <- species_range
    # Maximum aggregation since we don't want to loose any cells
    species_range50km <- raster::aggregate(species_range50km, fact = 5, na.rm = TRUE, fun = raster::modal)
    species_range50km <- alignRasters(species_range50km,globalgrid50km, method = "ngb",func = raster::modal, cl = FALSE)
    
    # Assert that some data is there
    assert_that(
      raster::cellStats(species_range,'max') == 1,
      raster::cellStats(species_range50km,'max') == 1
    )

    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction)
    # Clamp and convert to Integer
    species_range <- raster::clamp(species_range, lower = 0, upper = 1)
    species_range <- round( species_range * 1000 )
    
    # 10km
    # Save outputs  
    writeGeoTiff(species_range,
                 fname = out_10km,
                 dt = "INT2S"
    ) 
    
    # Now for 50km as well
    # Remove fractional cover
    species_range50km <- (species_range50km * globalgrid_fraction50km)
    # Clamp and convert to Integer
    species_range50km <- raster::clamp(species_range50km, lower = 0, upper = 1)
    species_range50km <- round( species_range50km * 1000 )
    
    # 50km
    writeGeoTiff(species_range50km,
                 fname = out_50km,
                 dt = "INT2S"
    ) 
    rm(species_range,species_range50km)
  } else {return(NULL)}
}
saveRDS(o, paste0("/media/martin/data/features/","newPlantPPM.rds") )
stopifnot( all( sapply(o, is.null) ) )

# ------------------------------------------------------- #
#### New Point data ####
# Get all range estimates for all species represented.
ll <- ll_points

# New outpath
output_path_10km = paste0("/media/martin/data/features/newPlantPoints_","10km")
output_path_50km = paste0("/media/martin/data/features/newPlantPoints_","50km")
if(!dir.exists(output_path_10km)) dir.create(output_path_10km)
if(!dir.exists(output_path_50km)) dir.create(output_path_50km)

# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = str_split(basename(ll),"__",simplify = T)[,2], fname = as.character(ll) ,stringsAsFactors = FALSE) 
plants$binomial <- str_replace(plants$binomial,"_points","")
plants$binomial <- tools::file_path_sans_ext(plants$binomial)

# Detect subspecies and correct subspecies to species
plants$subspecies <- Reduce("|", list(str_detect(plants$binomial,'subsp'), str_detect(plants$binomial,'var.')) )
plants_subspecies <- plants %>% dplyr::filter(subspecies == TRUE)
plants_subspecies$binomial <- paste(str_split(plants_subspecies$binomial,"_",simplify = T)[,1],str_split(plants_subspecies$binomial,"_",simplify = T)[,2])
# Join back
plants <- bind_rows( 
  plants %>% dplyr::filter(subspecies == FALSE),
  plants_subspecies
)
rm(plants_subspecies)

if(process_only_powo_species){
  plants <- matchName2POWO(plants, powo_checklist, powo_checklist_acc)
  # Also limit powo dist to the target ids
  powo_dist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/dist_NatureMap.txt') %>% 
    dplyr::filter(db_id %in% powo_checklist_acc$db_id ) %>% 
    dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}
data("tdwg")

# ------------------------ #
o <- foreach::foreach( row = 1:nrow(plants),
                       .inorder = F,
                       .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                       .export = c("ll","globalgrid","plants_unique",
                                   "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                                   "output_path_10km","output_path_50km","powo_dist","tdwg",
                                   "alignRasters","writeGeoTiff")
) %dopar% {
  out_10km = paste0(output_path_10km,"/",plants$binomial[row],".tif")
  out_50km = paste0(output_path_50km,"/",plants$binomial[row],".tif")
  if((!file.exists(out_10km)) & (!file.exists(out_50km)) ){
    species_range <- raster(plants$fname[row])
    
    myLog("Preparing species ",plants$binomial[row])
    
    # Reproject 
    species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
    
    # Alternative approach. Turn the grid cell into a point
    df <- as.data.frame(species_range,xy=TRUE)
    df <- df[,c('x','y',names(species_range))]
    df <- df[which(df[,3] == 1),]
    if( nrow(df)==0) {return(NULL)} # Empty raster
    
    # Make point    
    coordinates(df) <- ~x+y
    proj4string(df) <- proj4string(globalgrid)
    
    # Rasterize to the target grid
    species_range <- raster::rasterize(df,globalgrid,field = names(df), fun = 'first')
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid)
    )
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == plants$db_id[row])
    
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction)
    species_range <- round( species_range * 1000 )
    
    # 10km
    # Save outputs  
    writeGeoTiff(species_range,
                 fname = out_10km,
                 dt = "INT2S"
    ) 
    # --- #
    # Now for 50km
    # Rasterize to the target grid
    species_range <- raster::rasterize(df,globalgrid50km,field = names(df), fun = 'first')
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid50km)
    )
    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction50km)
    species_range <- round( species_range * 1000 )
    
    # 50km
    # Save outputs  
    writeGeoTiff(species_range,
                 fname = out_50km,
                 dt = "INT2S"
    ) 
    rm(species_range,df)
  } else { return(NULL) }
}
saveRDS(o, paste0("/media/martin/data/features/","NewPlantPoints.rds") )
stopifnot( all( sapply(o, is.null) ) )
gc()
raster::removeTmpFiles(.5)

# ------------------------------------------------------- #
#### New Bagging data ####
# Get all range estimates for all species represented in the NW
ll <- ll_rangebag
assert_that(length(ll)>0)

# New outpath
output_path_10km = paste0("/media/martin/data/features/newPlantRangebags_","10km")
output_path_50km = paste0("/media/martin/data/features/newPlantRangebags_","50km")
if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }

# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = str_split(basename(ll),"__",simplify = T)[,2],
                     type = str_split(basename(ll),"__",simplify = T)[,1],
                     fname = as.character(ll) ,stringsAsFactors = FALSE) %>% 
  # Filter to specific type
  dplyr::filter(type == "X0.165")

# Detect subspecies and correct subspecies to species
plants$subspecies <- Reduce("|", list(str_detect(plants$binomial,'subsp'), str_detect(plants$binomial,'var.')) )
plants_subspecies <- plants %>% dplyr::filter(subspecies == TRUE)
plants_subspecies$binomial <- paste(str_split(plants_subspecies$binomial,"_",simplify = T)[,1],str_split(plants_subspecies$binomial,"_",simplify = T)[,2])
# Join back
plants <- bind_rows( 
  plants %>% dplyr::filter(subspecies == FALSE),
  plants_subspecies
)
rm(plants_subspecies)

if(process_only_powo_species){
  plants <- matchName2POWO(plants, powo_checklist, powo_checklist_acc)
  # Also limit powo dist to the target ids
  powo_dist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/dist_NatureMap.txt') %>% 
    dplyr::filter(db_id %in% powo_checklist_acc$db_id ) %>% 
    dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}

# And the number of duplicated species
species_duplicates <- names(which(table(plants$binomial)>1)) 
# Number of duplicates
myLog("Number of globally duplicated species: ",length( species_duplicates ) )

# ------------------------ #
o <- foreach::foreach( row = 1:nrow(plants),
                       .inorder = F,
                       .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                       .export = c("ll","globalgrid","plants",
                                   "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                                   "output_path_10km","output_path_50km","powo_dist","tdwg",
                                   "alignRasters","writeGeoTiff")
) %dopar% {
  out_10km = paste0(output_path_10km,"/",plants$binomial[row],".tif")
  out_50km = paste0(output_path_50km,"/",plants$binomial[row],".tif")
  if((!file.exists(out_10km)) & (!file.exists(out_50km)) ){
    species_range <- raster(plants$fname[row])
    
    myLog("Preparing species ",plants$binomial[row])
    
    # Reproject
    species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
    
    # Alternative approach. Turn the grid cell into a point dataset. This generally seems to be faster than via raster
    df <- as.data.frame(species_range,xy=TRUE)
    df <- df[,c('x','y',names(species_range))]
    df <- df[which(df[,3] == 1),]
    if( nrow(df)==0) {return(NULL)} # Empty raster
    
    # Make point    
    coordinates(df) <- ~x+y
    proj4string(df) <- proj4string(globalgrid)
    
    # Rasterize to the target grid
    species_range <- raster::rasterize(df,globalgrid,field = names(df), fun = 'first')
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid)
    )
    
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == plants$db_id[row])
    
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    
    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction)
    species_range <- round( species_range * 1000 )
    
    # 10km
    # Save outputs  
    writeGeoTiff(species_range,
                 fname = out_10km,
                 dt = "INT2S"
    ) 
    # --- #
    # Now for 50km
    # Rasterize to the target grid
    species_range <- raster::rasterize(df,globalgrid50km,field = names(df), fun = 'first')
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid50km)
    )
    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction50km)
    species_range <- round( species_range * 1000 )
    
    # 50km
    # Save outputs  
    writeGeoTiff(species_range,
                 fname = out_50km,
                 dt = "INT2S"
    ) 
    rm(species_range,df)
  } else { return(NULL) }
}
saveRDS(o, paste0("/media/martin/data/features/","NewPlantRangebags.rds") )
assert_that(all( sapply(o, is.null) ))

stopCluster(cl);stopImplicitCluster()
stop("All done")