run_on_cluster = TRUE
# Purpose of this script is to create a copy of every range we got
# from which artifical land-surface types are masked out
library(tidyverse)
library(raster)
library(rgdal)
library(assertthat)
library(NatureMapRangeModelling)
library(gdalUtils)
library(fst)
library(fasterize)
library(exactextractr)
library(doParallel)
library(stringr)
source("src/000_ConvenienceFunctions.R")
withoutExtension <- function(x) { tools::file_path_sans_ext( x ) }
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
if(run_on_cluster){
  assessPOWOcoverage <- function(range, powo_extract, tdwg, verbose = TRUE) {
    # Import those function as needed
    pkg <- c('assertthat','raster')
    if(!all(pkg %in% installed.packages())) {
      stop(paste0("Package(s) ", pkg[which(!pkg %in% installed.packages())] ," are necessary for this function to work."))
    } else{
      require(assertthat);require(raster)
    }
    # Security checks
    stopifnot(
      assert_that(
        is.data.frame(powo_extract),
        has_name(powo_extract,"db_id"),
        length( unique(powo_extract$db_id)  ) == 1,
        class(range) == "RasterLayer",
        nlayers(range) == 1,
        class(tdwg) == "SpatialPolygonsDataFrame",
        has_name(tdwg,"LEVEL2_COD")
      )
    )
    
    # Subset and make a mask for the TDWG regions
    if(verbose) { print('-- Creating POWO mask --') }
    tdwg_sub <- subset(tdwg,LEVEL2_COD %in% unique(powo_extract$region_code_l2))
    tdwg_sub <- spTransform(tdwg_sub,CRSobj = CRS(proj4string(range)) )
    tdwg_sub <- rasterize(tdwg_sub,range,field = 'LEVEL2_COD')
    
    stopifnot(
      assert_that(
        compareRaster(tdwg_sub,range)
      )
    )
    # Now inverse mask out the labeled data
    if(verbose) { print('-- Masking out areas outside the TDWG --') }
    ras_zon <- suppressWarnings( raster::zonal(range, tdwg_sub, fun = 'max', na.rm = TRUE) )
    
    if( verbose) { print('-- Test coverage and remove otherwise --') }
    if( nrow(ras_zon) == 0 ) { return(NULL) }
    if( all(ras_zon[,2]>0) ){
      # All zones covered by at least one patch
      return(range)
    } else {
      return(NULL)
    }
    # Clean up
    rm(ras_zon,tdwg_sub)
    return(range)
  }
} else { assessPOWOcoverage <- NatureMapRangeModelling::assessPOWOcoverage }

# ------------------------------------- #
#### [Artifical masking] - Masking the plant data with the artifical mask ####
# ------------------------------------- #
# Plan:
# Remove from all vertebrates without habitat affil. the unsuitable habitat
# Assess coverage for plant species
cores <- max(c( parallel::detectCores() - 1, 6))
if(run_on_cluster){
  input_path = "/h/u046/jung/features" # Input path
  output_path = "/h/u046/jung/features_esh_artifical" # Output path
  art_mask = paste0( '/h/u046/jung/Artifical_lvl1_','10km','_mollweide_11082020.tif')
  art_mask50km = paste0( '/h/u046/jung/Artifical_lvl1_','50km','_mollweide_11082020.tif')
  powo_path = '/h/u046/jung/WCSP_KEW/'
} else{
  input_path = "/media/martin/data/features" # Input path
  output_path = "/media/martin/data/features_esh_artifical" # Output path
  art_mask = paste0( '/media/martin/data/Artifical_lvl1_','10km','_mollweide_11082020.tif') # Updated here
  art_mask50km = paste0( '/media/martin/data/Artifical_lvl1_','50km','_mollweide_11082020.tif') # Updated here     
  powo_path = '/media/martin/data/raw/WCSP_KEW/'
}
if(!dir.exists(output_path)) {dir.create(output_path)}
assert_that( file.exists(art_mask),
             dir.exists(input_path),
             dir.exists(output_path),
             dir.exists(powo_path)
             )

# Load POWO
powo_checklist <- data.table::fread(paste0(powo_path,'checklist_species.txt')) %>% 
  dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
powo_checklist_acc <- powo_checklist %>% dplyr::filter(taxon_status_description == 'Accepted')
powo_dist <- data.table::fread(paste0(powo_path,'dist_NatureMap.txt')) %>% 
  dplyr::filter(db_id %in% powo_checklist_acc$db_id )

# Grep folders with BIEN plant data
ld <- list.dirs(input_path,full.names = T)
ld <- ld[which(str_detect(ld, '10km'))] # Only target resolution

# Determine which folders have plant species in them from the name
plant_folders <- ld[union( which(str_detect(ld,'Plant')), which(str_detect(ld,'plant') ) )] # Check that all plant folders are detected
plant_folders <- c(plant_folders, ld[which(str_detect(ld,'PLANT')) ])

# Load mask 
artifical <- raster(art_mask)
artifical50km <- raster(art_mask50km)
assert_that(cellStats(artifical,"max")<=1000 )

# Reference grids
ref10km <- raster('data/globalgrid_mollweide_10km.tif')
ref50km <- raster('data/globalgrid_mollweide_50km.tif')

# Check and compare against reference grid
artifical <- projectRaster(artifical, crs = crs(ref10km))
artifical <- alignRasters(artifical, ref10km,method = 'bilinear',func = 'mean',cl = FALSE)
artifical50km <- projectRaster(artifical50km, crs = crs(ref50km))
artifical50km <- alignRasters(artifical50km, ref50km,method = 'bilinear',func = 'mean',cl = FALSE)

rm(ref10km,ref50km)

# Fire up cluster for parallel processing
if(file.exists('log_art.txt')){ file.remove('log_art.txt')}
cl <- parallel::makeCluster(cores, type = 'FORK',logfile = 'log_art.txt')
registerDoParallel(cl)
data("tdwg")  

# Loop through all folders and mask out in parallel
for( folder in ld ){
  # Create output folders if not already existing
  output_path_10km <- paste0(output_path,"/", str_remove(basename( folder ),"_10km"),"_10km" )
  output_path_50km <- paste0(output_path,"/", str_remove(basename( folder ),"_10km"),"_50km" )
  if(!dir.exists(output_path_10km)) { dir.create(output_path_10km) }
  if(!dir.exists(output_path_50km)) { dir.create(output_path_50km) }
  
  # List files
  lf <- list.files( folder, "*.tif",full.names = TRUE )
  assert_that(length(lf) > 0)
  
  if(folder %in% plant_folders){
    # Plant folder detected
    df <- data.frame(binomial = str_replace_all(tools::file_path_sans_ext(basename(lf))," ","_"),
                       fname = as.character(lf) ,stringsAsFactors = FALSE) 
    df$binomial <- str_replace(df$binomial,"\\.","_")
    df$binomial <- str_replace_all(df$binomial, " ", "_")
    df <- matchName2POWO(df, powo_checklist, powo_checklist_acc)
    pf <- TRUE # Check where powo coverage should be a criterion
  } else { 
    # Build data.frame
    df <- data.frame(binomial = str_replace_all(tools::file_path_sans_ext(basename(lf))," ","_"),
                     fname = lf, db_id = NA, stringsAsFactors = FALSE)
    pf <- FALSE
  }
  
  myLog('Processing ',basename(folder))
  # Subtract the fractional cover layer from each BIEN occurence map in parallel
  o <- foreach(row = 1:nrow(df),
               .inorder = FALSE,
               .errorhandling = 'pass',
               .export = c("output_path_10km","output_path_50km","writeGeoTiff","artifical","artifical50km","pf",'powo_dist','df','assessPOWOcoverage',
                           'tdwg', "alignRasters"),
               .packages=c('raster','assertthat','stringr','rgdal')
  ) %dopar% {
    out_10km = paste0(output_path_10km,"/",basename(df$fname[row]))
    out_50km = paste0(output_path_50km,"/",basename(df$fname[row]))
    if(file.exists(out_10km) & file.exists(out_50km)){ return(NULL) }
    
    print( df$binomial[row] )
    
    ras <- raster::raster(df$fname[row])
    
    if(!compareRaster(ras,  artifical,stopiffalse = F)){
      ras <- alignRasters(ras,artifical,method = "bilinear",func = raster::mean, cl = F)
    }
    # Sec check
    assert_that( compareRaster( ras,  artifical) )
    
    # Assess that there is enough coverage if those are plant species
    if(pf){
      # Get POWO subset
      powo_extract <- powo_dist %>% dplyr::filter(db_id == df$db_id[row])
      if(nrow(powo_extract)==0){ return(NULL)} # No distribution data, ignore
      # Check Coverage
      ras2 <- assessPOWOcoverage(range = ras,powo_extract = powo_extract,tdwg = tdwg)
      if(is.null(ras2)) { return(NULL)} else { ras <- ras2;rm(ras2)}
    }
    
    ras.copy = ras # Security copy
    # Remove patches at 10km
    ras = (ras - artifical)
    ras <- raster::clamp(ras, lower = 0, upper = 1000) # Ensure both lower and upper values are removed
    ras[ras == 0] <- NA
    # Write output
    if(cellStats(ras,'max')==Inf || cellStats(ras,'max')==-Inf){
      # Masking removed all habitat. Use range instead to retain the species
      ras = ras.copy
    }
    writeGeoTiff(ras, out_10km,'INT2U') # Save o  utputs and clear 
    
    # Now aggregate to 50km
    ras50km <- ras
    # Maximum aggregation since we don't want to loose any cells
    ras50km <- raster::aggregate(ras50km, fact = 5, na.rm = TRUE, fun = mean)
    ras50km <- alignRasters(ras50km,artifical50km, method = "bilinear",func = mean, cl = FALSE)
    ras50km[ras50km == 0] <- NA
    # Write output 
    writeGeoTiff(ras50km, out_50km,'INT2U') 
    # Done
    return(NULL)
  }
  saveRDS(o, paste0(output_path,"/", str_remove(basename( folder ),"_10km"),".rds") ) # Save log
  
} # Folder Loop done
stopCluster(cl);stopImplicitCluster()
stop('Done!')
