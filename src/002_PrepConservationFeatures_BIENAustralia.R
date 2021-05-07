# ------------------------------------------------------- #
# The preparation script for the BIEN Australia data  
# Creates rasterized species range layers per species
# Using a template grid as target

# 11/2019: Update
# Now also use POWO to (a) consider only target species in POWO as features and (b) remove unconnected patches in 
# non-native regions
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

# How many accepted species in POWO don't have distribution data?
table(powo_checklist_acc$hasdistribution)/ sum(table(powo_checklist_acc$hasdistribution))

# BIEN path
bien_path <- paste0("/media/martin/data/raw/Plantdata_BIENAustralia/")

# Output path
output_path_10km = paste0("/media/martin/data/features/BIENAustraliaplants_","10km")
output_path_50km = paste0("/media/martin/data/features/BIENAustraliaplants_","50km")
if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }

# Some error checks
stopifnot(
  assert_that(
    file.exists(bien_path),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif")),
    dir.exists(output_path_10km),
    dir.exists(output_path_50km)
  )
)

# ------------------------------------------------------- #
# Get all range estimates for all species represented.
ll <- c(
  list.files(paste0(bien_path,"AUSTRALIANbinary_tp05rangemaps/"),"*.tif",recursive = T,full.names = TRUE)
)
stopifnot( assert_that(length(ll)>0) )

# Make a data.frame of all inputs and split basenames
plants <- data.frame(binomial = tools::file_path_sans_ext(basename(ll)), fname = as.character(ll), stringsAsFactors = FALSE) 
plants$binomial <- str_replace(plants$binomial,"\\.","_")
plants$binomial <- str_replace_all(plants$binomial," ","_")
  
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
myLog("Number of duplicated species: ",length( species_duplicates ) )

# ------------------------ #
# Make a global land area map out of the grid
globalgrid <- raster(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif"))
globalgrid <- setMinMax(globalgrid)
globalgrid[globalgrid>=1] <- 1
globalgrid50km <- raster(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif"))
# And the fraction cover
globalgrid_fraction <- raster(paste0(gridpath,"globalgrid_",projection,"_","10km",".tif")) / 1000
globalgrid_fraction50km <- raster(paste0(gridpath,"globalgrid_",projection,"_","50km",".tif")) / 1000
# ------------------------ #
# Subset to POWO
if(process_only_powo_species){
  plants <- matchName2POWO(plants, powo_checklist, powo_checklist_acc)
  # Also limit powo dist to the target ids
  powo_dist <- powo_dist %>% dplyr::filter(db_id %in% unique(plants$db_id))
  myLog( n_distinct(plants$binomial), ' species found in POWO')
}

# Process first the non-duplicated species!
plants_unique <- plants %>% dplyr::filter(binomial %notin% species_duplicates)

data("tdwg")
# Fire up cores
file.remove('log.txt')
cl <- makePSOCKcluster(cores,outfile="log.txt")
registerDoParallel(cl = cl)

# Start the process
o <- foreach::foreach( row = 1:nrow( plants_unique ),
                       .errorhandling = 'pass',
                  .packages = c("fasterize","raster","stringr","velox","NatureMapRangeModelling",
                                "assertthat","gdalUtils","landscapemetrics"),
                  .export = c("ll","globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                              "output_path_10km","output_path_50km","powo_dist","tdwg",
                              "plants_unique","alignRasters","writeGeoTiff")
) %dopar% {
  out_10km = paste0(output_path_10km,"/",plants_unique$binomial[row],".tif")
  out_50km = paste0(output_path_50km,"/",plants_unique$binomial[row],".tif")
  if((!file.exists(out_10km)) & (!file.exists(out_50km)) ){
    species_range <- raster(plants_unique$fname[row])
    
    myLog("Preparing species: ",plants_unique$binomial[row])
    
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == plants_unique$db_id[row])
    
    # Reproject
    species_range <- projectRaster(species_range,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE )
    
    # Extend to resolution of global grid
    species_range <- raster::extend(species_range, globalgrid)
    
    # Align resolution
    species_range[species_range==0] <- NA
    species_range <- alignRasters(species_range,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
    
    # Remove unconnected non-native patches
    # ----- #
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    if(raster::cellStats(species_range,'max') == 0) { return(NULL)} # Filtering removed all habitat patches (unlikely)
    # ----- #
    
    # Also aggregate to 50km using a maximum 
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
  }
}
# Save output list to check later
saveRDS(o, paste0("/media/martin/data/features/","BIENAustraliaplants.rds") )
stopifnot( all( sapply(o, is.null) ) )

# Then the duplicated species
plants_nonunique <- plants %>% dplyr::filter(binomial %in% species_duplicates)

myLog("Now processing the duplicated plant entries. Total: ", n_distinct(plants_nonunique$binomial))

# New outpath
output_path_10km = paste0("/media/martin/data/features/BIENAustraliaplants_duplicates_","10km")
output_path_50km = paste0("/media/martin/data/features/BIENAustraliaplants_duplicates_","50km")
if(!dir.exists(output_path_10km)) dir.create(output_path_10km)
if(!dir.exists(output_path_50km)) dir.create(output_path_50km)

o <- foreach::foreach( row = 1:nrow(plants_nonunique),
                  .inorder = T,
                  .packages = c("fasterize","raster","stringr","velox","assertthat","gdalUtils","landscapemetrics","NatureMapRangeModelling"),
                  .export = c("ll","globalgrid","globalgrid50km","globalgrid_fraction","globalgrid_fraction50km",
                              "output_path_10km","output_path_50km","powo_dist","tdwg",
                              "plants_nonunique","alignRasters","writeGeoTiff")
) %dopar% {
  out_10km = paste0(output_path_10km,"/",plants_nonunique$binomial[row],".tif")
  out_50km = paste0(output_path_50km,"/",plants_nonunique$binomial[row],".tif")
  if((!file.exists(out_10km)) & (!file.exists(out_50km)) ){
    myLog("Preparing species ",plants_nonunique$binomial[row])
  
    # Get all species in the list that have the species name and Load them as stack
    ras_names <- plants$fname[which( plants$binomial == plants_nonunique$binomial[row] )]  
    ras_names <- ras_names[grep("xml",basename(ras_names),invert = T)] # Remove XML files if existing
  
    species_range <- raster::stack()
    for(sr in ras_names){
      sr <- raster::raster(sr)
      
      # Reproject
      sr <- projectRaster(sr,crs = CRS( proj4string(globalgrid) ),method = 'ngb',alignOnly = FALSE  )
      
      # Extend to resolution of global grid
      sr <- raster::extend(sr, globalgrid)
      
      # Align resolution
      sr[sr==0] <- NA
      sr <- alignRasters(sr,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
      species_range <- addLayer(species_range, sr)
    }
    # Aggregate both subsets
    species_range <- sum(species_range,na.rm = T)
    species_range[species_range>=1] <- 1
    species_range[species_range==0] <- NA
  
    # Assert that some data is there
    assert_that(
      raster::cellStats(species_range,'max') == 1
    )
    
    # Get POWO subset
    powo_extract <- powo_dist %>% dplyr::filter(db_id == plants_nonunique$db_id[row])
    
    # Align resolution
    species_range <- alignRasters(species_range,globalgrid, method = "ngb",func = raster::modal, cl = FALSE)
    
    # Remove unconnected non-native patches
    # ----- #
    if(nrow(powo_extract)>0){
      species_range <- NatureMapRangeModelling::removeNonNativePatches(
        range = species_range,tdwg = tdwg,powo_extract = powo_extract, verbose = FALSE
      )
    }
    if(raster::cellStats(species_range,'max') == 0) { return(NULL)} # Filtering removed all habitat patches (unlikely)
    # ----- #
    
    # Also aggregate to 50km using a maximum 
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
  } else { return(NULL) }
}
# Save output list to check later
saveRDS(o, paste0("/media/martin/data/features/","BIENAustraliaplants_duplicates.rds") )
stopifnot( all( sapply(o, is.null) ) )

stopCluster(cl);stopImplicitCluster()
stop("DONE!")
