# Purpose of this script is to load and aggregate the GEE processed 
# ESH files. This is done in R rather than python for convenience
# (not nessarily for speed)
library(tidyverse)
library(tidyr)
library(raster)
library(rgdal)
library(assertthat)
library(NatureMapRangeModelling)
library(gdalUtils)
library(stars)
library(fst)
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

what =  '10km' # Resolutions to process
what_alt = '11039m' # In new export called accordingly 
cores = parallel::detectCores() - 1 # For parallel processing

# Path parameters 
input_path  = paste0("/media/martin/data/temp_gee_aoh") # Changed to base path
output_path = "/media/martin/data/features_esh" # Output path
if(!dir.exists(output_path)) dir.create(output_path)
output_path_aoo = "/media/martin/data/features" # Output path to base AOO
species_path = paste0("data/speciesID_table_",what,".fst")
powo_path = '/media/martin/data/raw/WCSP_KEW/checklist_species.txt'
species_list_gard = "data/newGardREP_IDs.csv" # From Andy. Species from GARD that are not in IUCN
species_list_bgci = "data/PLANTS_NON_iucn_gdb.csv"
# FIXME: New Plant list from 28/08/2020
#full_plant_list = 'data/PLANTAE_IUCN_2019_2.csv'
marine_habs_only = "data/marine_habs_only_rl2019vw_200817.csv" # Marine habitats only

# Load global grids 
global_raster <- raster(paste0('data/globalgrid_','mollweide','_',what,'.tif'))
land_mask <- global_raster; land_mask[land_mask>0] <- 1
global_raster_wgs84 <- raster(paste0('data/globalgrid_','wgs84','_',what,'.tif'))

ignore_existing = FALSE # Should existing files be loaded and merged in? 
multiplier = 1000 # Multiplier! (Nov. 2019 = 1000000)
# Global 
temp_path     = "/tmp/mergesave"
if(!dir.exists(temp_path)) dir.create(temp_path)

# Output log file
log_file <- data.frame(iucn_id_no = character(), bandid = numeric(),fname = character(), comment = character(),stringsAsFactors = FALSE)
log_file_path <- paste0(input_path,'/','log_AOH_processing_',what,'.csv')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ------------------------------------- #
# Read in species files
# ------------------------------------- #
species_list  <- read.fst(species_path) %>% dplyr::filter(data %in% c(
  "IUCN Mammal", "IUCN Reptiles", "IUCN Amphibians", "IUCN Bird", "IUCN Plants", "Kew IUCN Plants", "Kew BGCI Plants", "GARD Shai" 
))
species_list_gard <- read_csv( species_list_gard ) %>% dplyr::filter(id_no_type == 'GARD') %>% 
  dplyr::select(id_no, binomial)
# Make a master species list
master_list1 <- species_list %>% dplyr::filter(data != "GARD Shai") %>% 
  dplyr::select(iucn_id_no, binomial, family, fname, data, folder)
master_list2 <- species_list %>% dplyr::filter(data == 'GARD Shai') %>% 
  dplyr::select(iucn_id_no, binomial, family, fname, data, folder)
# Remove reptiles already present in IUCN
master_list2 <- master_list2[-which( str_to_lower(master_list2$binomial) %in% str_to_lower(master_list1$binomial) ),]
# Assemble the master list
master_list <- bind_rows(master_list1,master_list2)
# Cleanup
rm(master_list1,master_list2, species_list_gard)

# And for the plants
species_list_bgci <- read_csv( species_list_bgci )
powo_checklist <- data.table::fread(powo_path) %>% dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
plant_ids <- species_list_bgci %>% dplyr::mutate(fname = paste0(str_replace_all(BINOMIAL," ","_"),".tif"),
                                                 data = 'GEE_PLANTAE', folder = paste0('PLANTAE_BGCI_',what) ) %>% 
  dplyr::rename(iucn_id_no = ID_NO, binomial = BINOMIAL) %>% 
  dplyr::select(iucn_id_no, binomial, fname, data,folder)

# Species with marine habitats only
# marine_habs_only <- read_csv(marine_habs_only)
# Won't currently exclude. Many of these coastal birds still requiere coastal areas

myLog( "Number of species processed ", length(  unique(master_list$iucn_id_no) ))
# Security checks
stopifnot(
  assert_that(
    anyNA(master_list$folder) == FALSE,
    anyNA(master_list$iucn_id_no) == FALSE,
    anyNA(master_list$fname) == FALSE,
    anyNA(master_list$binomial) == FALSE,
    dir.exists(input_path)
  )
)

# ------------------------------------- #
#### Process all files ####
# ------------------------------------- #
# Load, split and preprocess all file

# List of folders to process
ld <- list.dirs(input_path)
if(length(ld)>1) ld <- ld[-1] # Without base path
#ld <- ld[str_detect(ld,paste0("_",what))] # Get only those files of a given resolution

# Last checks
stopifnot(
  assert_that(
  length(ld) > 0,
  length( list.files(temp_path,'*.tif') ) == 0,  # Template folder is empty
  msg = 'Necessary input data missing or temporary folder is full of files'
  )
)

# Fire up cluster #
if(file.exists(paste0(output_path,'/log_ESH.txt'))){ file.remove(paste0(output_path,'/log_ESH.txt'))}
cl <- parallel::makeCluster(cores, type = 'FORK',log = paste0(output_path,'/log_ESH.txt'))
registerDoParallel(cl)

# ------------------------------------- #
# -------- Start of processing -------- #
# [Processing] 2020-08-29 08:53:24 | Mosaicing of stats_AOH_v1_multiband_habmap_v003_REPTILIA_imageScale_11039m_statsScale_11039m_global.csv done. Cleaning up now...
#Error: nrow(band_species) not equal to nb
# "stats_AOH_v1_multiband_habmap_v003_plantaeBGCI_imageScale_11039m_statsScale_11039m_global.csv"

# Now for each folder
for(folder in ld){
  # First load tifs
  list_of_tifs <- list.files(folder,pattern = "*.tif",full.names = T)
  list_of_tifs <- list_of_tifs[has_extension(list_of_tifs,"tif")]
  
  # Then band names
  list_of_csvs <- list.files(folder, pattern = "*.csv",full.names = TRUE)
  list_of_csvs <- list_of_csvs[has_extension(list_of_csvs,"csv")]
  
  # Sec check that the lists are not empty
  stopifnot( assert_that(length(list_of_tifs)>0, length(list_of_csvs)>0 ) )
  
  # Starting
  myLog('Currently processing ', basename(folder))
  
  for( csv in list_of_csvs ){

    # Load associated csv
    band_species <- read_csv( csv ) %>% tibble::rownames_to_column('bandid') %>% 
      # Separate id coloumn into groups
      tidyr::separate(id,into = c("TGrouping","Status","iucn_id_no","season"),sep = '_') %>% 
      dplyr::select(bandid, filter:season)
    stopifnot(n_distinct(band_species$filter)==1)
    
    # Which file to use. Use Andy's list for the plants for now
    if( length(grep('plantaeBGCI',basename(csv))) >0){ ids <- plant_ids} else {  ids <- master_list }
    
    # Join with masterlist to get the species names
    band_species <- left_join(band_species, ids %>% 
                                dplyr::select(iucn_id_no, binomial) %>% distinct() %>% 
                                mutate(iucn_id_no = as.character(iucn_id_no)), by = 'iucn_id_no')
    
    if(anyNA(band_species$binomial)){
      log_file <- bind_rows(log_file,
        band_species[ which(is.na(band_species$binomial)),] %>% dplyr::select(iucn_id_no,bandid) %>% 
          mutate(fname = basename(csv), bandid = as.numeric(bandid), comment = 'Not found in species master list.')
      )
    }
    stopifnot(anyDuplicated(band_species$bandid)==0) # Problem with joining to lost?
    
    # Construct output filenames
    band_species <- mutate(band_species, ofname = paste0( as.character(iucn_id_no),"__", str_replace_all(str_to_lower(binomial),' ',"_"),".tif" ) )

    # Get the associated tifs with the same Taxonomic group
    target_tifs <- list_of_tifs[
      str_detect(basename(list_of_tifs), unique(band_species$filter) )
    ]
    # Get number of bands in layer
    nb <- nlayers( raster::stack(target_tifs[1]) )
    stopifnot(assert_that(nrow(band_species)==nb, length(target_tifs) >0 ))

    # -- #
    # Do some other house keeping on folder metadata
    tgroup <- unique(band_species$filter)
    tthreat <- basename(folder) # Batch
    tres <- what 
  
    # Create target Directory if don't exist if not alredy existing
    out_dir = paste0( output_path,"/",paste0(tgroup,"_", tres) )
    if(!dir.exists( out_dir )) { dir.create( out_dir )}
    
    # Create a folder for GARD
    out_dir_gard = paste0( output_path,"/",paste0("reptiliaGARD","_", tres) )
    if(!dir.exists( out_dir_gard )) { dir.create( out_dir_gard )}

    # Check whether existing files should be ignored?
    if(ignore_existing) { 
      # Construct outputs and check if they are there
      if(unique(band_species$TGrouping) == 'reptiliaGARD') { out_dir_new <- out_dir_gard} else 
      {out_dir_new <- out_dir} # Set a new outdir if it is a GARD reptile
      list_of_fnames <- paste0(out_dir_new, "/", band_species$ofname)
      # Skip to next csv in loop
      if(all(file.exists(list_of_fnames)) == TRUE) {next()}
      }
    
    # Clean up warnings
    last.warning <- NULL
    # Security checks
    stopifnot(assert_that(
      tgroup != "", tres != "", tthreat != "",
      dir.exists(out_dir),
      length(list_of_tifs) > 0),
      nrow(band_species) == nb,
      any(!duplicated(band_species$bandid)) # Any duplicated band names
  )
  myLog("     --> Processing ", nb ," species")
  # For each file split them up
  for( fraction_file in target_tifs) {
    myLog(" --> Reading in file ", basename( fraction_file ))
    # Load as stars array
    x <- read_stars(fraction_file,proxy = FALSE)
    
    # Security check
    if(!is.null(stars::st_get_dimension_values(x,'band'))) { # If check if there is only one band
      stopifnot( assert_that(length( stars::st_get_dimension_values(x,'band') ) == nb) )
    }
    
    # Split all files in raster
    for(i in 1:nb){
      if(band_species$bandid[i] == 'constant'){ next()}
      # Construct temporary file name | Link back to species_id |
      # Single file in folder
      fbname <- paste0( temp_path,"/", tgroup,
                        "__",band_species$bandid[i],"__", # band id is critical. This is how they are merged
                        tools::file_path_sans_ext(basename(fraction_file)),
                        ".tif" )
      # Slice up array
      if(nb>1){
        # Multiple bands
        y <- x %>% slice(index = i, along = "band")
      } else { y <- x } # Single band
      
      # Write output
        write_stars(y,
                    dsn = fbname,
                    options = c("COMPRESS=DEFLATE"),
                    NA_value = 0, # Set zero to NA to allow easier aggregation
                    type = ifelse( multiplier>65534, 'UInt32', 'UInt16')
       )
      }
    } # All bands processed
  rm(x,y, fbname) # Clean
  
  myLog("...Starting mosaicing and output creation...")
  
  # All the individual tiles for mosaicing
  target_files = list.files(temp_path,"*.tif",full.names = TRUE)
  target_files <- target_files[has_extension(target_files,"tif")]
  
  # Start the great parallel mosaicing
  o <- foreach(sn = iter( band_species$bandid ),
               .inorder = FALSE,
               .errorhandling = 'pass',
               .export = c("out_dir","temp_path","target_files","writeGeoTiff",'multiplier','out_dir_gard',
                           "land_mask","global_raster","global_raster_wgs84","alignRasters", 'withoutExtension','ignore_existing'),
               .packages = c('raster','assertthat','gdalUtils','stringr'),
               .verbose = FALSE) %dopar% {
                 if(sn == 'constant'){ return(NULL)}
                 # Get band $species by matching to master table
                 species_name <- band_species %>% dplyr::filter( bandid == sn )
                 if(is.na(species_name$binomial)){ return(NULL) } # Species not found in master file, saved earlier in log. skip
                 # output name
                 if(nrow(species_name)==0) {
                   # Security copy / use id instead for filename
                   fname <- paste0(out_dir, "/", sn,".tif")
                   fname.vrt <- paste0(temp_path, "/", sn,".vrt")
                 } else {
                   if(nrow(species_name) == 1 && unique(species_name$TGrouping) == 'reptiliaGARD') { out_dir_new <- out_dir_gard} else 
                     {out_dir_new <- out_dir} # Set a new outdir if it is a GARD reptile
                   fname <- paste0(out_dir_new, "/", species_name$ofname)
                   fname.vrt <- paste0(temp_path, "/", withoutExtension( species_name$ofname ),".vrt")
                 }
                 
                 # Check whether existing files should be ignored?
                 if(ignore_existing && file.exists(fname)) { return(NULL) }
                 
                 # Get all files to mosaic
                 list_of_tifsm <- target_files[which(str_detect(target_files,paste0('__',sn,'__')))]
                 if(length(list_of_tifsm)>0){
                   # Create temporary vrt
                   gdalbuildvrt(
                     gdalfile = list_of_tifsm,
                     output.vrt = fname.vrt,
                     resolution = "highest",
                     vrtnodata = -9999,
                     r = "average",
                     separate = FALSE,
                     ot = ifelse( multiplier>65534, 'UInt32', 'UInt16'),
                     overwrite = TRUE
                   )
                   ras <- raster::raster(fname.vrt)
                   # Extend to extend of global grid
                   ras <- raster::extend(ras, global_raster_wgs84)
                   # Crop to global grid (get's rid of Antarctica and other bits)
                   ras <- raster::crop(ras, global_raster_wgs84, snap = 'near')
                   # Reproject to mollweide
                   ras <- projectRaster(ras,crs = proj4string(global_raster),method = 'bilinear',alignOnly = FALSE)
                   # Replace zeros with NA
                   ras[ras == 0] <- NA
                   # And align because grid shape differences
                   ras <- alignRasters(ras,global_raster,method = 'bilinear',func = mean, cl = FALSE)
                   
                   # Load output file to memory if it already exists and average
                   if(file.exists(fname)){ 
                     ras.existing <- raster(fname) 
                     ras.existing[ras.existing<1] <- NA # overwrite NA to be safe
                     ras <- mean(stack(ras, ras.existing),na.rm = TRUE)
                     ras[ras==0] <- NA
                     rm(ras.existing)
                   }
                   # Mask with land area mask
                   ras <- raster::mask(ras, land_mask)
                   
                   # Convert / round to nearest Integer
                   ras <- round(ras)
                   
                   # Write output again
                   writeGeoTiff(ras, fname,'INT2U') # Save outputs and clear 
                   rm(ras) # Delete output
                   file.remove(fname.vrt) # Delete vrt file
                   return(NULL)
                 } # End of stacking loop
               } # End of foreach loop
  # Check failed tasks
  if(!all(sapply(o, is.null)) ){
    # Some tasks failed, save information on why in the log file
    failed_ids <- which(!sapply(o, is.null))
    # o[failed_ids,]
    log_file <- bind_rows(
      log_file,
      band_species %>% dplyr::filter(bandid %in% failed_ids) %>% dplyr::select(iucn_id_no, bandid) %>% 
        dplyr::mutate(bandid = as.numeric(bandid), fname = basename(csv), comment = as.character(o[failed_ids]))
    )
    rm(failed_ids)
  }

  gc() # Garbarge collecting
  print('\n')
  myLog("Mosaicing of ", basename(csv), " done. Cleaning up now...")
  
  # Clean up temp folder
  unlink( target_files,force = TRUE)
  } # End of csv file loop
} # End of folder loop
                  
# Now query missing species and attach to the log
ll_csvs <- list.files(input_path, '*.csv',recursive = TRUE, full.names = TRUE) %>% 
  purrr::map_df(., .f = function(x){
    read_csv( x ) %>% tibble::rownames_to_column('bandid') %>% 
      # Separate id coloumn into groups
      tidyr::separate(id,into = c("TGrouping","Status","iucn_id_no","season"),sep = '_') %>% 
      dplyr::select(bandid, filter:season) %>% mutate(csvname = x)
  })

# Join with masterlist to get the species names
ll_csvs <- left_join(ll_csvs, master_list %>% 
                            dplyr::select(iucn_id_no, binomial) %>% distinct() %>% 
                            mutate(iucn_id_no = as.character(iucn_id_no)), by = 'iucn_id_no')

if(anyNA(ll_csvs$binomial)){
  log_file <- bind_rows(log_file,
                        ll_csvs[ which(is.na(ll_csvs$binomial)),] %>% dplyr::select(iucn_id_no,bandid,csvname) %>% 
                          mutate(fname = basename(csvname), bandid = as.numeric(bandid), comment = 'Not found in species master list.')
  )
}
# Drop missing species
ll_csvs <- ll_csvs %>% drop_na(binomial)

# Construct output filenames
ll_csvs <- mutate(ll_csvs, ofname = paste0( as.character(iucn_id_no),"__", str_replace_all(str_to_lower(binomial),' ',"_"),".tif" ) )

# Then load all tif files and check which ones are there
ll_tifs <- list.files(output_path, '*.tif',recursive = TRUE, full.names = TRUE)
ll_tifs <- ll_tifs[has_extension(ll_tifs,'tif')]

# Missing species are those that were not written to the output folder
log_file <- bind_rows(log_file,
                      ll_csvs[which( ll_csvs$ofname %notin% basename(ll_tifs)),] %>% dplyr::select(iucn_id_no,bandid,csvname) %>% 
                        mutate(fname = basename(csvname), bandid = as.numeric(bandid), comment = 'Empty band(s). No output written.')
)

# Write log file
write_csv(log_file,log_file_path)
stopImplicitCluster()
stop('!DONE!')

# --------------------------------------- #
# Some species got not matched by name to the previous list
# Requested an additional ID list from WCMC to rename the existing ones
additional_species_list_bgci = read.csv("data/plants_without_name_IUCN_GDB.csv",stringsAsFactors = FALSE) # Joined in by Corinna to missing names
# Make a filename similar to IUCN plants
additional_species_list_bgci$new.fname <- paste0(additional_species_list_bgci$binomial_1,".tif")
additional_species_list_bgci$full_fname <- str_replace( additional_species_list_bgci$full_fname, paste0("PLANTAE_", "10km"), paste0("PLANTAE_", what))
additional_species_list_bgci$folder <- paste0("PLANTAE_", what)
assert_that( all( file.exists(additional_species_list_bgci$full_fname) ) )

# Rename all those files
file.rename(from = additional_species_list_bgci$full_fname,
            to = paste0( dirname(additional_species_list_bgci$full_fname), "/",additional_species_list_bgci$new.fname )
            )
rm(additional_species_list_bgci)
# --- #
# Check for remaening fit with IUCN species list
iucn_list <- read_fst(paste0('data/speciesID_table_',what,'.fst')) # The global biodiversity path lists
# All species with ID
list_of_tifs <- list.files(output_path,recursive = T,full.names = T)
list_of_tifs <- list_of_tifs[which(str_detect(list_of_tifs,'id_'))] # Get only those are named id
miss <- data.frame(full_fname = list_of_tifs, wrong_fname = basename(list_of_tifs),stringsAsFactors = F) %>% 
  dplyr::mutate(iucn_id_no = as.numeric( str_remove(tools::file_path_sans_ext(basename(list_of_tifs)),'id_'))  ) %>% 
  dplyr::left_join(., iucn_list %>% dplyr::select(iucn_id_no, fname ) %>% distinct() ) %>% 
  dplyr::mutate(folder = dirname(list_of_tifs))

file.rename(from = miss$full_fname,
            to = paste0( miss$folder,'/',miss$fname )
)
rm(miss,list_of_tifs,iucn_list)

# --------------------------------------- #
#### Aggregate 10km to 50km ####
# --------------------------------------- #
# Take the existing grid and aggregate to reference grid
# No mean aggregation. Ensure that suitable habitat is <= total area in PU

# Reference grid
globalgrid_50km <- raster('data/globalgrid_mollweide_50km.tif')

# Target folders
ld <- list.dirs(output_path)
ld <- ld[str_detect(ld,'10km')] # Only the 10km directiories

if(file.exists('log.txt')){ file.remove('log.txt')}
cl <- parallel::makeCluster(cores, type = 'FORK',logfile = 'log.txt')
registerDoParallel(cl)

# Loop over each folder
for(folder in ld){
  print(basename(folder))
  
  # New folder name and create
  odname <- paste0(
    dirname(folder),'/',
    str_remove(basename(folder),'_10km'), '_50km'
    )
  dir.create(odname,showWarnings = FALSE)
  
  tifs <- list.files(folder,'*.tif',full.names = TRUE)
  tifs <- tifs[has_extension(tifs,'tif')]
  stopifnot( assert_that(length(tifs)>0) )
  
  # Create data.frame of input tifs, outputfile
  df <- data.frame(fname = tifs, 
                   iucn_id_no = str_split(basename(tifs),"_",simplify = T)[,1],
                   ofname = paste0(odname, '/', basename(tifs)),
                   stringsAsFactors = FALSE)
  # Link file names with available EOO masks
  df <- dplyr::left_join(df, master_list %>% dplyr::select(iucn_id_no, folder, fname) %>% dplyr::rename(inname = fname) %>% 
                           dplyr::mutate(iucn_id_no = as.character(iucn_id_no)), by = 'iucn_id_no' )
  #stopifnot( assert_that(all(is.na(df$inname)==FALSE)) )
  #df <- df %>% mutate(aoo_mask = paste0("/media/martin/data/features/",df$folder,"/",df$inname) )

  # Start processing
  o <- foreach(row = 1:nrow(df),
               .inorder = FALSE,
               .errorhandling = 'pass',
               .export = c("writeGeoTiff","df","globalgrid_50km",
                           "alignRasters", "master_list",'withoutExtension'),
               .packages=c('raster','assertthat','stringr','rgdal'),
               .verbose = TRUE) %dopar% {
                 sub <- df[row,]
                 myLog('Processing ', basename(sub$fname) )
                 # Skip if already existing
                 if(file.exists(sub$ofname)) {return(NULL)}
                 
                 # Load AOH and EOO raster
                 aoh <- raster::raster( sub$fname )
                 # if(!is.na(sub$aoo_mask) || !file.exists(sub$aoo_mask)){
                 #   eoo <- raster::raster( sub$aoo_mask )
                 #   eoo.mask <- eoo > 0 # Make mask
                 #   if(compareRaster(aoh,eoo) == FALSE){
                 #     aoh <- alignRasters(aoh, eoo,method = 'bilinear',func = mean, cl = FALSE)
                 #   }
                 #   # Mask AOH with eoo.mask
                 #   aoh <- raster::mask(aoh, eoo.mask)
                 #   # Overwrite aoh layer
                 #   # writeGeoTiff(aoh, sub$full.fname,'INT2U')
                 # }
                 
                 # Overwrite missing data with NA to be sure
                 aoh[aoh < 1] <- NA
                 
                 # Now aggregate the layer    
                 aoh_50km <- alignRasters(aoh, globalgrid_50km, method = 'bilinear',func = mean, cl =FALSE)

                 # Security check that data is still present
                 if(cellStats(aoh_50km,'max')==-Inf) { return(NULL) } # Correct later then
                 # Write output again
                 writeGeoTiff(aoh_50km, sub$ofname, 'INT2U')
                 rm(aoh,sub,aoh_50km)
                 return(NULL)
               } # End of masking foreach loop
  # Security checks
  o[which(!sapply(o, is.null))]
  
}

ld <- list.dirs(output_path)
ld <- ld[str_detect(ld,'10km')] # Only the 10km directiories
files_10km <- list.files(path = ld,pattern = "*.tif",full.names = TRUE)
files_10km <- files_10km[has_extension(files_10km,'tif')]
ld <- list.dirs(output_path)
ld <- ld[str_detect(ld,'50km')] # Only the 10km directiories
files_50km <- list.files(path = ld,pattern = "*.tif",full.names = TRUE)
files_50km <- files_50km[has_extension(files_50km,'tif')]
# Differentiate
x <- setdiff2(basename(files_10km), basename(files_50km))

# Copy for investigation
dir.create('/home/martin/EmptyRasters_10km/',showWarnings = FALSE)
file.copy(files_10km[which(basename(files_10km) %in% x$X_not_in_Y)],
          paste0('/home/martin/EmptyRasters_10km/',
                 basename(dirname(files_10km[which(basename(files_10km) %in% x$X_not_in_Y)])),
                 basename(files_10km[which(basename(files_10km) %in% x$X_not_in_Y)])))
file.remove(files_10km[which(basename(files_10km) %in% x$X_not_in_Y)])

# --------------------------------------- #
#### Masking of extracted AOH with EOO ####
# --------------------------------------- #
# Plan:
# Loop through all AOH files and mask them with the EOO
# If file not existing, insert EOO instead
# Save EOO estimate in separate folder to ease processing
# Important to keep track of files that got masked and files that are (still) missing

# Target resolution 
what = '10km'

# List of AOH folders
ld <- list.dirs(output_path);ld <- ld[-1]
ld <- ld[which(str_detect(ld,what))] # Filter to target resolution
assert_that(length(ld) > 0, length(list.files(ld)) > 0,
            exists('master_list')) # Security checks

# Get all files and match them to a respective range
list_of_tifs <- list.files(ld,pattern = "*.tif",full.names = TRUE,recursive = TRUE)

# Construct matching data.frame for gee extracts  
gee_extracts <-  data.frame(fname = basename(list_of_tifs), full.fname = list_of_tifs, folder = basename(dirname(list_of_tifs)), what = what,
                            stringsAsFactors = FALSE ) %>% 
  dplyr::mutate(binomial = str_to_lower(tools::file_path_sans_ext(fname) ))
gee_extracts$binomial <- ifelse( str_count(gee_extracts$binomial,"_") == 2,
                                 paste0( str_split(gee_extracts$binomial,"_",simplify = T)[,2],"_",str_split(gee_extracts$binomial,"_",simplify = T)[,3]),
                                 gee_extracts$binomial)
gee_extracts$binomial <- str_to_sentence(str_replace(gee_extracts$binomial,"_"," "))
assert_that(anyNA(gee_extracts$binomial)==FALSE,anyNA(gee_extracts$folder) == FALSE)
# -- #
# Now match in the species names first by filename
EOO_found <- inner_join(gee_extracts %>% dplyr::select(-binomial), master_list, by = c('fname','folder') ) %>% 
  dplyr::mutate(eoo.folder = folder, eoo.fname = fname)
EOO_notfound <- anti_join(gee_extracts, master_list, by = c('fname','folder')) # Basically all plants
# Try and match the notfound species by binomial name only
EOO_found2 <- inner_join(EOO_notfound %>% dplyr::filter(folder == paste0('PLANTAE_',what) ),
                         master_list %>% dplyr::filter(folder == paste0('PLANTAE_',what) ) %>% 
                           dplyr::rename(eoo.folder = folder) %>% 
                           distinct() %>% dplyr::rename(eoo.fname = fname),
                         by = 'binomial')
# And try for different feature names
EOO_found3 <- inner_join(EOO_notfound %>% dplyr::filter(folder == paste0('PLANTAE_BGCI_',what) ),
                         master_list %>% dplyr::filter(folder %in% c( paste0('KEW_PLANTAE_BGCI_',what), paste0('KEW_PLANTAE_IUCN_',what) )) %>% 
                           dplyr::rename(eoo.folder = folder) %>% 
                           distinct() %>% dplyr::rename(eoo.fname = fname),
                         by = 'binomial')
# Generally by binomial
EOO_found4 <- inner_join(EOO_notfound %>% dplyr::filter(fname %notin% EOO_found2$fname,fname %notin% EOO_found3$fname),
                         master_list %>% dplyr::rename(eoo.folder = folder) %>% 
                           distinct() %>% dplyr::rename(eoo.fname = fname),
                         by = 'binomial')
# Combine all those names that could be matched
gee_extracts2 <- bind_rows(
    EOO_found,EOO_found2,EOO_found3,EOO_found4
  )
# Those are still missing.
x <- gee_extracts[which(gee_extracts$fname %notin% gee_extracts2$fname),]
gee_extracts <- gee_extracts2

rm(EOO_found,EOO_found2,EOO_found3,EOO_found4,EOO_notfound,x,gee_extracts2)

# Rename those files here on disk
# FIXME: Would have been better if that were corrected during the EOO step.
x = file.exists(paste0(output_path_aoo,"/",gee_extracts$eoo.folder,"/",gee_extracts$eoo.fname))
gee_extracts[!x,'eoo.fname'] <- str_replace_all(gee_extracts[!x,'eoo.fname'],' ','_')
rm(x)

# Make sure that all AOH & EOO files exist and everything is correct
assert_that(
  anyNA(gee_extracts$fname) == FALSE, anyNA(gee_extracts$full.fname) == FALSE,
  anyNA(gee_extracts$eoo.fname) == FALSE,   anyNA(gee_extracts$eoo.folder) == FALSE,
  all( file.exists(gee_extracts$full.fname) ),
  all( file.exists( paste0(output_path_aoo,"/",gee_extracts$eoo.folder,"/",gee_extracts$eoo.fname)) )
)

# Fire up cluster for parallel processing
if(file.exists('log_ESH.txt')){ file.remove('log_ESH.txt')}
cl <- parallel::makeCluster(cores, type = 'FORK',logfile = 'log_ESH.txt')
registerDoParallel(cl)

# Start processing
o <- foreach(row = 1:nrow(gee_extracts),
               .inorder = FALSE,
               .errorhandling = 'pass',
               .export = c("writeGeoTiff","gee_extracts","output_path_aoo",
                           "alignRasters", "master_list",'withoutExtension'),
               .packages=c('raster','assertthat','stringr','rgdal'),
               .verbose = TRUE) %dopar% {
                 
                 sub <- gee_extracts[row,]
                 myLog('Processing ', sub$binomial, ' at ', sub$what, ' resolution' )
                 
                 # Load AOH and EOO raster
                 aoh <- raster::raster( sub$full.fname )
                 eoo <- raster::raster(  paste0(output_path_aoo,"/",sub$eoo.folder,"/",sub$eoo.fname) ) # Load EOO  (must be present!)
                 eoo.mask <- eoo > 0
                 # Security check
                 assert_that(
                   cellStats(eoo.mask,'max') == 1,
                   compareRaster(aoh,eoo)
                 )
                 
                 # Mask AOH with eoo.mask
                 aoh <- raster::mask(aoh, eoo.mask)
                 # Overwrite aoh layer
                 writeGeoTiff(aoh, sub$full.fname,'INT2U')
                 
                 # Next ensure that AOH is never greater than EOO
                 aoh <- min(aoh,eoo)
                 # Security check that data is still present
                 if(cellStats(aoh,'max')==-Inf) { return(NULL) } # Correct later then
                 # Write output again
                 writeGeoTiff(aoh, sub$full.fname,'INT2U')
                 return(NULL)
               } # End of masking foreach loop
  # Security checks
  assert_that(
    all(sapply(o, is.null))
    )

stopCluster(cl);
stopImplicitCluster()
stop('DONE! Finished masking')

# ------------------------------------- #
#### Correct non-native and sufficient coverage of plants only ####
# This is done for the artificially masked species but not yet for the AOH data
warning('Ensure that a backup exist. Script below alters files and removes them!')
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

# Load POWO
powo_checklist <- data.table::fread(paste0("/media/martin/data/raw/WCSP_KEW/",'checklist_species.txt')) %>% 
  dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
powo_checklist_acc <- powo_checklist %>% dplyr::filter(taxon_status_description == 'Accepted')
powo_dist <- data.table::fread(paste0("/media/martin/data/raw/WCSP_KEW/",'dist_NatureMap.txt')) %>% 
  dplyr::filter(db_id %in% powo_checklist_acc$db_id )

# Get PLANT folders and files
plant_folders <- list.dirs(output_path); plant_folders <- plant_folders[str_detect(plant_folders,"PLANTAE")]
plant_files <- list.files(plant_folders,recursive = TRUE, full.names = TRUE)

# Check and match names to POWO
plant_df <- data.frame(binomial = str_replace_all(tools::file_path_sans_ext(basename(plant_files))," ","_"),
                 fname = as.character(plant_files) ,stringsAsFactors = FALSE) 
plant_df$binomial <- str_replace(plant_df$binomial,"\\.","_")
plant_df$binomial <- str_replace_all(plant_df$binomial, " ", "_")
plant_df <- matchName2POWO(plant_df, powo_checklist, powo_checklist_acc)

# --- #
# Warning: The code below deletes all the files not found in POWO!
warning('MAKE A BACKUP!')
p2r = plant_files[which(plant_files %notin% plant_df$fname)]
sapply(p2r, file.remove)
# --- #

# Fire up cluster #
if(file.exists('log_ESH.txt')){ file.remove('log_ESH.txt')}
cl <- parallel::makeCluster(cores, type = 'FORK',log = 'log_ESH.txt')
registerDoParallel(cl)

# Check POWO coverage and delete if not enough
o <- foreach(row = 1:nrow(plant_df),
             .inorder = FALSE,
             .errorhandling = 'pass',
             .export = c("plant_df","plant_folders","tdwg"),
             .packages=c('raster','assertthat','stringr','rgdal')
) %dopar% {
  print( plant_df$binomial[row] )
  
  ras <- raster::raster(plant_df$fname[row])
  
  if(cellStats(ras,'max')==0 || cellStats(ras,'max')==-Inf){ 
    rm(ras)
    file.remove(plant_df$fname[row])
    return(NULL)
    }
  # Assess that there is enough coverage if those are plant species
  # Get POWO subset
  powo_extract <- powo_dist %>% dplyr::filter(db_id == plant_df$db_id[row])
  if(nrow(powo_extract)==0){ return(NULL)} # No distribution data, ignore
  # Check Coverage
  ras2 <- assessPOWOcoverage(range = ras,powo_extract = powo_extract,tdwg = tdwg)
  if(is.null(ras2)) { rm(ras,ras2);  file.remove(plant_df$fname[row]) } 
  return(NULL)
}

# ------------------------------------- #
####    Gdal built and aggregate     ####
# ------------------------------------- #

# First mosaic and warp all images together

folder = list_of_tifs[2]
setwd( paste0(input_path,"/",folder) )
gdalbuildvrt(
  gdalfile = list.files(".","*.tif"),
  output.vrt = paste0(folder,".vrt"),
  resolution = "highest",
  srcnodata	= NULL,
  vrtnodata = -9999,
  separate = FALSE,
  ot = 'Byte',
  r = 'nearest',
  overwrite = TRUE
)
# Now warp  
gdalwarp(
  srcfile = paste0(folder,".vrt"),
  dstfile = paste0(folder,".tif"),
  r = 'near',
  srcnodata = -9999,
  dstnodata	= -9999,
  ot = 'Byte',
  multi = TRUE,
  co = c("NUM_THREADS=ALL_CPUS","COMPRESS=DEFLATE"),
  verbose=TRUE
)
