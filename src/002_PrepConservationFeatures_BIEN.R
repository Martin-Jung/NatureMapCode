# ------------------------------------------------------- #
# The preparation script for the BIEN data  
# Creates rasterized species range layers per species
# Using a template grid as reference

# 11/2019: Update
# Now also use POWO to (a) consider only target species in POWO as features and (b) remove unconnected patches in 
# non-native regions. Furthermore harmonize folder structure and file names for plants and aggregate all data together if duplicates occur.
library(raster)
library(fasterize)
library(rgdal)
library(gdalUtils)
library(NatureMapRangeModelling) # Customized package not available as part of this release
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
  plants <- dplyr::left_join(plants, powo_checklist_acc %>% dplyr::select(acc_full_name_without_family,db_id) %>% distinct() %>% 
                               dplyr::mutate(acc_full_name_without_family = str_to_lower(str_replace_all(acc_full_name_without_family," ","_"))),
                             by = c("binomial" = "acc_full_name_without_family") )
  
  # Now also match the synonms
  # Now check with all the synonym's
  plants_syn <- dplyr::filter(plants,is.na(db_id)) %>% dplyr::select(-db_id) %>% 
    dplyr::left_join(., powo_checklist %>% dplyr::select(full_name_without_family,acc_full_name_without_family) %>% 
                       dplyr::mutate(acc_full_name_without_family = str_to_lower( str_replace_all(acc_full_name_without_family," ","_"))),
                     by = c("binomial" = "full_name_without_family")) %>% 
    drop_na(acc_full_name_without_family)  # Drop all remaining without accepted species name
  # Join in db_id of accepted species again
  plants_syn <- dplyr::left_join( plants_syn, powo_checklist_acc %>% dplyr::select(acc_full_name_without_family,db_id) %>% distinct() %>% 
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

# Path to BIEN modelled plant ranges
dir_bien_ppm <- c(
  "/media/martin/data/raw/Plantdata_BIENAustralia/AUSTRALIANbinary_tp05rangemaps/",
  "/media/martin/data/raw/BIEN/bien41_outputs/PlantsUpdated72619/NWPlants/_outputsAllMerged/PPM/BinaryMaps/",
  "/media/martin/data/raw/BIEN/bien41_outputs/zipped_CEIBA_plant_AT/_outputs/PPM/BinaryMaps/",
  "/media/martin/data/raw/BIEN/bien41_outputs/zipped_CEIBA_plant_IM/_outputs/PPM/BinaryMaps/"
  )
dir_bien_points <- c(
  "/media/martin/data/raw/BIEN/bien41_outputs/PlantsUpdated72619/NWPlants/_outputsAllMerged/Points/BinaryMaps/"
)
dir_bien_rangebag <- c(
  "/media/martin/data/raw/BIEN/bien41_outputs/PlantsUpdated72619/NWPlants/_outputsAllMerged/RangeBag/BinaryMaps/"
)
assert_that(
  all( sapply(dir_bien_ppm, dir.exists) ),
  all( sapply(dir_bien_points, dir.exists) ),
  all( sapply(dir_bien_rangebag, dir.exists) )
)
# POWO data
process_only_powo_species <- TRUE # Only process species that can be assigned to POWO taxonomy
powo_checklist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/checklist_species.txt') %>% 
  dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
powo_checklist_acc <- powo_checklist %>% dplyr::filter(taxon_status_description == 'Accepted')
powo_dist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/dist_NatureMap.txt') %>% 
  dplyr::filter(db_id %in% powo_checklist_acc$db_id )

# Some error checks
stopifnot(
  assert_that(
    all( dir.exists( dir_bien_ppm ) ),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif"))
  )
)

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
#### BIEN PPM models ####
# Get all range estimates for all species represented.
ll <- list.files(dir_bien_ppm,recursive = TRUE, full.names = TRUE)

# Output path
output_path_10km = paste0("/media/martin/data/features/BIEN_PLANTAE_PPM_","10km")
output_path_50km = paste0("/media/martin/data/features/BIEN_PLANTAE_PPM_","50km")
if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }

# Some error checks
stopifnot(
  assert_that(
    length(ll)>0,
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif")),
    dir.exists(output_path_10km),
    dir.exists(output_path_50km)
  )
)

# ------------------------------------------------------- #
# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = ifelse(str_detect(basename(tools::file_path_sans_ext(basename(ll))),"__"),
                                       str_split(basename(tools::file_path_sans_ext(basename(ll))),"__",simplify = T)[,2],
                                       tools::file_path_sans_ext(basename(ll))
                                       ), fname = as.character(ll) ,stringsAsFactors = FALSE)
plants$binomial <- str_replace(plants$binomial,"\\.","_")
plants$binomial <- str_replace(plants$binomial," ","_")
stopifnot( 
  assert_that( length(which(plants$binomial == '')) == 0,
               anyNA(plants$binomial) == FALSE
               )
)
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
  plants$binomial <- str_to_lower(plants$binomial)
  plants_withid <- matchName2POWO(plants, powo_checklist, powo_checklist_acc) %>% dplyr::select(binomial, db_id)
  stopifnot( assert_that( all( plants_withid$binomial %in% plants$binomial ) ) )
  plants <- left_join(plants_withid, plants, by = 'binomial') # This now also has the duplicated species in here
  # Also limit powo dist to the target ids
  powo_dist <- powo_dist %>% dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}
data("tdwg")

# Fire up cores
if(file.exists('BIEN_log.txt')) {file.remove('BIEN_log.txt')}
cl <- makeCluster(cores)#, outfile = "BIEN_log.txt")
registerDoParallel(cores = cores)

# And the number of duplicated species
species_duplicates <- names(which(table(plants$binomial)>1)) 
# Number of duplicates
myLog("Number of globally duplicated species: ",length( species_duplicates ) )

# Process
# ------------------------ #
o <- foreach::foreach( sn = iter( plants$binomial ),
                    .inorder = F,
                    .errorhandling = 'pass',
                    .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                    .export = c("ll","globalgrid","plants",
                                "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                                "output_path_10km","output_path_50km","powo_dist","tdwg",
                                "alignRasters","writeGeoTiff")
  ) %dopar% {
    # Subset of plant species
    sub <- subset(plants, binomial == sn)
    out_10km = paste0(output_path_10km,"/",sn,".tif")
    out_50km = paste0(output_path_50km,"/",sn,".tif")
    if(file.exists(out_10km) & file.exists(out_50km)){ return(NULL) } # Both output files already existing
    
    if(nrow(sub) >1) {
      # Get all species in the list that have the species name and Load them as stack
      ras_names <- unique(sub$fname)
      ras_names <- ras_names[assertthat::has_extension(ras_names,'tif')] # Only tif files
      
      species_range <- raster::stack()
      for(fn in ras_names){
        sr <- raster::raster(fn)
        
        # Reproject
        sr <- projectRaster(sr,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
        
        # Extend to resolution of global grid
        sr <- raster::extend(sr, globalgrid)
        
        # Align resolution
        sr[sr==0] <- NA
        sr <- alignRasters(sr,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
        species_range <- addLayer(species_range, sr)
        rm(sr)
      }
      # Aggregate both subsets
      species_range <- sum(species_range,na.rm = T)
      species_range[species_range>=1] <- 1
      species_range[species_range==0] <- NA
      
      # Assert that some data is there
      assert_that(
        compareRaster(species_range, globalgrid),
        raster::cellStats(species_range,'max') == 1
      )
    } else {
      species_range <- raster(sub$fname)
    }
    if(is.infinite(cellStats(species_range,'max'))) { return(NULL)} # 33005 is empty?
    
    print(paste0("Preparing species = ",sn) )
    
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == sub$db_id[1])
    
    # Reproject
    species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE)
    
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
    return(NULL)
}
stopifnot( all( sapply(o, is.null) ) )
saveRDS(o, paste0("/media/martin/data/features/","BIENplants.rds") )
gc()
raster::removeTmpFiles(.5)
# ------------------------------------------------------- #
#### BIEN Point data ####
# Get all range estimates for all species represented.
ll <- list.files(dir_bien_points,recursive = TRUE, full.names = TRUE)

# New outpath paths
output_path_10km = paste0("/media/martin/data/features/BIEN_PLANTAE_POINTS_","10km")
output_path_50km = paste0("/media/martin/data/features/BIEN_PLANTAE_POINTS_","50km")
if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }

# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = ifelse(str_detect(basename(tools::file_path_sans_ext(basename(ll))),"__"),
                                       str_split(basename(tools::file_path_sans_ext(basename(ll))),"__",simplify = T)[,1],
                                       tools::file_path_sans_ext(basename(ll))
), fname = as.character(ll) ,stringsAsFactors = FALSE)
plants$binomial <- str_remove_all(plants$binomial, "__points")
plants$binomial <- str_remove_all(plants$binomial, "_points")
plants$binomial <- str_replace(plants$binomial,"\\.","_")
plants$binomial <- str_replace(plants$binomial," ","_")
plants$binomial <- str_to_lower(plants$binomial)

stopifnot( 
  assert_that( length(which(plants$binomial == '')) == 0,
               length( which(str_detect(plants$binomial,"points")) ) == 0,
               anyNA(plants$binomial) == FALSE
  )
)
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
  plants$binomial <- str_to_lower(plants$binomial)
  plants_withid <- matchName2POWO(plants, powo_checklist, powo_checklist_acc) %>% dplyr::select(binomial, db_id)
  stopifnot( assert_that( all( plants_withid$binomial %in% plants$binomial ) ) )
  plants <- left_join(plants_withid, plants, by = 'binomial') # This now also has the duplicated species in here
  # Also limit powo dist to the target ids
  powo_dist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/dist_NatureMap.txt') %>% 
    dplyr::filter(db_id %in% powo_checklist_acc$db_id ) %>% 
    dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}
data("tdwg")
# And the number of duplicated species (although this is irrelevant for points processing)
species_duplicates <- names(which(table(plants$binomial)>1)) 
# Number of duplicates
myLog("Number of globally duplicated species: ",length( species_duplicates ) )

# ------------------------ #
o <- foreach::foreach( sn = iter(plants$binomial),
                  .inorder = F,
                  .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                  .export = c("ll","globalgrid","plants",
                              "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                              "output_path_10km","output_path_50km","powo_dist","tdwg",
                              "alignRasters","writeGeoTiff")
                  ) %dopar% {
    sub <- subset(plants, binomial == sn)
    
    out_10km = paste0(output_path_10km,"/",sn,".tif")
    out_50km = paste0(output_path_50km,"/",sn,".tif")
    if(file.exists(out_10km) & file.exists(out_50km)){ return(NULL) }
    myLog("Preparing species ",sn)
    # Load in all the point rasters
    out <- data.frame()
    for(i in 1:nrow(sub)){
      species_range <- raster(sub$fname[i])
      # Reproject 
      species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
      # Alternative approach. Turn the grid cell into a point
      df <- as.data.frame(species_range,xy=TRUE)
      df <- df[,c('x','y',names(species_range))]
      df <- df[which(df[,3] == 1),]
      if( nrow(df)==0) {next()} else { # Empty raster
        out <- rbind(out, df)
      }
    }
    df <- out
    if(nrow(df)==0){ return(NULL)}
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
    powo_extract <- powo_dist %>% dplyr::filter(db_id == sub$db_id[1])
    
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    if(cellStats(species_range,'max') == -Inf){ return(NULL) }
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
}
saveRDS(o, paste0("/media/martin/data/features/","BIENplantspoints.rds") )
stopifnot( all( sapply(o, is.null) ) )
gc()
raster::removeTmpFiles(.5)

# ------------------------------------------------------- #
#### BIEN Range bagging data ####
# Get all range estimates for all species represented in the NW
ll <- list.files(dir_bien_rangebag,recursive = TRUE, full.names = TRUE)

# New outpath paths
output_path_10km = paste0("/media/martin/data/features/BIEN_PLANTAE_RANGEBAG_","10km")
output_path_50km = paste0("/media/martin/data/features/BIEN_PLANTAE_RANGEBAG_","50km")
if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }

# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = str_split(basename(ll),"__",simplify = T)[,2],
                     type = str_split(basename(ll),"__",simplify = T)[,1],
                     fname = as.character(ll) ,stringsAsFactors = FALSE) %>% 
  # Filter to specific type
  dplyr::filter(type == "X0.165") %>% dplyr::select(-type)
plants$binomial <- str_replace(plants$binomial,"\\.","_")
plants$binomial <- str_replace(plants$binomial," ","_")
plants$binomial <- str_to_lower(plants$binomial)

stopifnot( 
  assert_that( length(which(plants$binomial == '')) == 0,
               length( which(str_detect(plants$binomial,"X0")) ) == 0,
               anyNA(plants$binomial) == FALSE
  )
)
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
  plants$binomial <- str_to_lower(plants$binomial)
  plants_withid <- matchName2POWO(plants, powo_checklist, powo_checklist_acc) %>% dplyr::select(binomial, db_id)
  stopifnot( assert_that( all( plants_withid$binomial %in% plants$binomial ) ) )
  plants <- left_join(plants_withid, plants, by = 'binomial') # This now also has the duplicated species in here
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
o <- foreach::foreach( sn = iter( plants$binomial ),
                       .inorder = F,
                       .errorhandling = 'pass',
                       .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","NatureMapRangeModelling","landscapemetrics"),
                       .export = c("ll","globalgrid","plants",
                                   "globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                                   "output_path_10km","output_path_50km","powo_dist","tdwg",
                                   "alignRasters","writeGeoTiff")
) %dopar% {
  # Subset of plant species
  sub <- subset(plants, binomial == sn)
  out_10km = paste0(output_path_10km,"/",sn,".tif")
  out_50km = paste0(output_path_50km,"/",sn,".tif")
  if(file.exists(out_10km) & file.exists(out_50km)){ return(NULL) } # Both output files already existing
  
  if(nrow(sub) >1) {
    # Get all species in the list that have the species name and Load them as stack
    ras_names <- unique(sub$fname)
    ras_names <- ras_names[assertthat::has_extension(ras_names,'tif')] # Only tif files
    
    species_range <- raster::stack()
    for(fn in ras_names){
      sr <- raster::raster(fn)
      
      # Reproject
      sr <- projectRaster(sr,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
      
      # Extend to resolution of global grid
      sr <- raster::extend(sr, globalgrid)
      
      # Align resolution
      sr[sr==0] <- NA
      sr <- alignRasters(sr,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
      species_range <- addLayer(species_range, sr)
      rm(sr)
    }
    # Aggregate both subsets
    species_range <- sum(species_range,na.rm = T)
    species_range[species_range>=1] <- 1
    species_range[species_range==0] <- NA
    
    # Assert that some data is there
    assert_that(
      compareRaster(species_range, globalgrid),
      raster::cellStats(species_range,'max') == 1
    )
  } else {
    species_range <- raster(sub$fname)
  }
  if(is.infinite(cellStats(species_range,'max'))) { return(NULL)} # 33005 is empty?
  
  print(paste0("Preparing species = ",sn) )
  
  # Get POWO subset
  powo_extract <- powo_dist %>% dplyr::filter(db_id == sub$db_id[1])
  
  # Reproject
  species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE)
  
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
  return(NULL)
}
saveRDS(o, paste0("/media/martin/data/features/","BIENplantsrangebags.rds") )
stopifnot( all( sapply(o, is.null) ) )
gc()
raster::removeTmpFiles(.5)
stopCluster(cl);stopImplicitCluster()
stop("Points all done")