# Idea of this script is to create a global PU Table per PU id and Biodiversity feature (FID)
# Runs in parallel for a target resolution
if(dir.exists('/h/u046/jung/naturemap')) { setwd('/h/u046/jung/naturemap'); run_on_cluster = TRUE } else { run_on_cluster = FALSE }
# Number of cores and cluster_type for parallel processing 
cores <- max( parallel::detectCores() - 1, 7)
cluster_type <- "FORK" 
    
# Paths
tmp_dir <- tempdir() # set directory to store temporary files
main_path <- ifelse(run_on_cluster, '/h/u046/jung/','/media/martin/data/')

output_path <- "data" # Output path for database
fst_path <- paste0(main_path,'results') # Save locally for now
dir.create(fst_path,showWarnings = FALSE)
output_filename <- "table_puspecies" # The final output table

# Target resolution and projection
target_resolution = "50km"
target_projection = "mollweide"
target_range      = c("","_esh")[2]

biodiversity_path <- paste0( main_path, 'features', target_range ) # Where is the feature data stored
aoh_artifical_path <- paste0( main_path, 'features',target_range,'_artifical' ) # Where is the feature data stored
eoo_path <- paste0( main_path, 'features')
also_use_eoo <- TRUE # use the EOO in case no AOH is available
constrain_path <- paste0(main_path,'constraints/')
constrain_type <- c("modal","fraction")[2]

# Create a log file for parallel processing
logfile <- paste0("log_putable_",target_resolution,".txt")
if(file.exists(logfile)) file.remove(logfile) # Remove if already present

# ------------------------------------------------------ #
library(raster)
library(gdalUtils)
library(prioritizr)
library(tidyverse)
library(assertthat)
library(fst) 
library(doParallel)
source("src/000_ConvenienceFunctions.R")

# Input pu data
pu <- raster(paste0(output_path,"/globalgrid_",target_projection,"_",target_resolution,".tif"))
pu <- setMinMax(pu)

# Make a copy of the PU layer and insert the cell number as PU id
pu_id_raster <- pu # Make a copy to count up pu id's
pu_cells <- raster::Which(!is.na(pu), cells = TRUE)
pu_id_raster[pu_cells] <- seq_along(pu_cells)
# create planning unit data in tabular format
pu_data <- data.frame(id = seq_along(pu_cells),
                      cost = pu[pu_cells])
pu[pu_cells] <- 1 # Equal cost everywhere
assert_that(
  cellStats(pu,'max') == 1,
  all( between(pu_data$cost, 0, 1000) ),
  max(pu_data$id) < .Machine$integer.max
)

# Species #
# -------------------------------- # 
# Matchup table existing
if(file.exists(paste0(output_path,"/speciesID_table_",target_resolution,".fst"))){
  species_ids <- read.fst(paste0(output_path,"/speciesID_table_",target_resolution,".fst"))
} else{
  myLog("Matchup table not found. Need to create a new one")
  stop('No species matchup table found.')
  species_ids <- read.fst(paste0(output_path,"/speciesID_table_",target_resolution,".fst"))
}
stopifnot(
  assert_that( nrow(species_ids) > 0,
               class(pu) == "RasterLayer",
               dir.exists(biodiversity_path), dir.exists(fst_path))
)

# Duplicated Reptile species of different datasets
x1 <- species_ids[which(species_ids$data == "GARD Shai"),] # Get all reptiles from GARD not in IUCN reptiles
x2 <- species_ids[which(species_ids$data == "IUCN Reptiles"),] 
species_ids <- species_ids %>% dplyr::filter(iucn_id_no %notin% x1[which(x1$binomial %in% x2$binomial),'iucn_id_no']) # Remove them from the species_id list
rm(x1,x2)

# Slightly more complicated for the plants
# All plant data
plant_ids <- species_ids %>% dplyr::filter(data %in% c("IUCN Plants",
                                                       "Kew IUCN Plants", "Kew BGCI Plants",
                                                       "New Plants - PPM", "New Plants - Points", "New Plants - Rangebag",             
                                                       "BIEN Plantae PPM", "BIEN Plantae Points" ) 
                                           ) %>% dplyr::select(data,binomial, iucn_id_no)
myLog('Loaded ', nrow(plant_ids), ' plant species')

# Process by data source and in order of quality
plant_ids1_1 <- plant_ids %>% dplyr::filter(data %in% c("IUCN Plants"))
plant_ids1_2 <- plant_ids %>% dplyr::filter(data %in% c("Kew IUCN Plants")) %>% 
  dplyr::filter(binomial %notin% plant_ids1_1$binomial )
plant_ids1_3 <- plant_ids %>% dplyr::filter(data %in% c("Kew BGCI Plants")) %>% 
  dplyr::filter(binomial %notin% plant_ids1_1$binomial,binomial %notin% plant_ids1_2$binomial )
plant_ids1 <- bind_rows(plant_ids1_1,plant_ids1_2,plant_ids1_3);rm(plant_ids1_1,plant_ids1_2,plant_ids1_3)

# Then the new global models of species that are not already in the expert range data
plant_ids2_1 <- plant_ids %>% dplyr::filter(data %in% c("New Plants - PPM")) %>%  dplyr::filter(binomial %notin% plant_ids1$binomial )
plant_ids2_2 <- plant_ids %>% dplyr::filter(data %in% c("New Plants - Rangebag")) %>%  dplyr::filter(binomial %notin% plant_ids1$binomial,
                                                                                                     binomial %notin% plant_ids2_1$binomial)
plant_ids2_3 <- plant_ids %>% dplyr::filter(data %in% c("New Plants - Points")) %>%  dplyr::filter(binomial %notin% plant_ids1$binomial,
                                                                                                   binomial %notin% plant_ids2_1$binomial,
                                                                                                   binomial %notin% plant_ids2_2$binomial)
# Assemble
plant_ids2 <- bind_rows(plant_ids2_1,plant_ids2_2,plant_ids2_3); rm(plant_ids2_1,plant_ids2_2,plant_ids2_3)

# Then everything else not already in the new species modelled
plant_ids3_1 <- plant_ids %>% dplyr::filter(data %in% c("BIEN Plantae PPM")) %>% dplyr::filter(binomial %notin% plant_ids1$binomial,binomial %notin% plant_ids2$binomial)
plant_ids3_2 <- plant_ids %>% dplyr::filter(data %in% c("BIEN Plantae Points")) %>% dplyr::filter(binomial %notin% plant_ids1$binomial,binomial %notin% plant_ids2$binomial) %>% 
  dplyr::filter(binomial %notin% plant_ids3_1$binomial)
# Assemble
plant_ids3 <- bind_rows(plant_ids3_1,plant_ids3_2); rm(plant_ids3_1,plant_ids3_2)

plant_ids <- bind_rows(plant_ids1,plant_ids2,plant_ids3) %>% distinct();rm(plant_ids1,plant_ids2,plant_ids3)
myLog('Cleaned to ', nrow(plant_ids), ' from those species')
stopifnot( assert_that(anyDuplicated(plant_ids$iucn_id_no) == 0) )

# Now get the final set of species for the extraction
species_ids <- bind_rows(
  species_ids %>% dplyr::filter(data %in% c("IUCN Mammal", "IUCN Reptiles", "IUCN Amphibians", "IUCN Bird","GARD Shai") ),
  left_join(plant_ids, species_ids, by = c('data','binomial','iucn_id_no') )
)
rm(plant_ids)
# Remove the remaining duplicated IDs for now (2 in total)
species_ids <- species_ids %>% dplyr::filter(iucn_id_no %notin% species_ids[which(duplicated(species_ids$iucn_id_no)),'iucn_id_no'] )

# Also remove dolphins, wales, seal & dugongs for sure
species_ids <- species_ids %>% dplyr::filter( str_to_lower(species_ids$family) %notin% c('delphinidae','balaenopteridae','dugongidae',
                                                                                         'balaenidae','cetotheriidae','balaenopteridae',
                                                                                          'eschrichtiidae','monodontidae','phocoenidae',
                                                                                         'physeteridae','kogiidae','platanistidae',
                                                                                         'iniidae','pontoporiidae','ziphiidae',
                                                                                         'phocidae','odobenidae','otariidae'))

stopifnot( assert_that(anyDuplicated(species_ids$iucn_id_no) == 0) )

# ---------------------------------------------------- #
# Get all target main features from the target range
feature_files <- list.files(biodiversity_path,recursive = TRUE, full.names = TRUE)
feature_files <- feature_files[has_extension(feature_files,'tif')]
myLog('Loaded ', length(feature_files), ' for processing ', target_range)
# --- #
# Make data frame
features <- data.frame(fname = basename(feature_files), full.fname = feature_files, folder = basename(dirname(feature_files)),stringsAsFactors = FALSE) %>% 
  # Get only target resolution
  dplyr::filter(str_detect(folder, target_resolution)) %>% 
  dplyr::mutate(type = ifelse(target_range == '','Range','AOH'))
if(target_range=='_esh'){
  features <- features %>% mutate(iucn_id_no = str_split(fname,"__",simplify = T)[,1] ) %>% 
    # Subset to only those filenames names in the species_id file
    dplyr::filter(iucn_id_no %in% species_ids$iucn_id_no)
} else {
  features <- features %>% 
  # Subset to only those filenames names in the species_id file
  dplyr::filter(str_to_lower(str_replace_all(fname,' ','_'))[,] %in% str_to_lower(str_replace_all(species_ids$fname,' ','_')) )
}

# Any duplicates?
if(anyDuplicated(features$iucn_id_no)>0){
  # Remove the duplicates 3 Reptiles with assessed by both IUCN and GARD
  features <- features[-which(duplicated(features$iucn_id_no)),]
}  
# Security check
assert_that(nrow(features)>0);rm(feature_files)
# --- #

# Now get the artifical AOH alternative
if(target_range != '') {
  # Load the files in batches
  dd <- list.dirs(aoh_artifical_path)
  dd <- dd[-1];  dd <- dd[str_detect(dd,target_resolution)]
  stopifnot(
    assert_that(length(dd) == n_distinct(species_ids$folder))
  )
  # Load the files in batches
  feature_art_files <- vector()
  for(folder in dd){
    myLog('Reading in = ', basename(folder))
    ff <- list.files(folder,full.names = TRUE)
    ff <- ff[has_extension(ff,'tif')] 
    feature_art_files <- c(feature_art_files, ff)
    rm(ff)
  }
  stopifnot(
    assert_that(
      length(feature_art_files) > 190000
    )
  )
  # --- #
  # Make data frame as well
  features_art <- data.frame(fname = basename(feature_art_files), full.fname = feature_art_files,
                             folder = basename(dirname(feature_art_files)),stringsAsFactors = FALSE) %>% 
    # Get only target resolution
    dplyr::filter(str_detect(folder, target_resolution)) %>% 
    # Subset to only those filenames names in the species_id file
    dplyr::filter(str_to_lower(str_replace_all(fname,' ','_')) %in% str_to_lower(str_replace_all(species_ids$fname,' ','_')) ) %>% 
    dplyr::mutate(type = 'Artifical AOH')
  # Join in the id from species id file
  features_art <- dplyr::inner_join(
    species_ids %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))) %>% 
      dplyr::select(iucn_id_no, binomial,join.fname, folder, data) %>% distinct(),
    features_art %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))),
    by = c('folder','join.fname')
  ) %>% mutate(iucn_id_no = as.character(iucn_id_no)) %>% dplyr::select(-join.fname,-binomial,-data)
  
  # Security check
  stopifnot( assert_that(nrow(features_art)>0, n_distinct(features_art$folder) == n_distinct(species_ids$folder))  )
  rm(feature_art_files)
  # Make copy for later
  features_art_backup <- features_art
  # --- #
  # Now filter to only those artifial files that are not already in features
  features_art <- features_art %>% 
    dplyr::filter(iucn_id_no %notin% features$iucn_id_no )
  # Now append to the features list
  features <- bind_rows(
    features, features_art
  )
  rm(features_art)
  
  # Print a summary
  #table(features$type,features$folder)
} else {
  features_art_backup <- c()
}

# Join in the id from species id file
if(!has_name(features,'iucn_id_no')){
  features <- dplyr::inner_join(
    species_ids %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))) %>% 
      dplyr::select(iucn_id_no, binomial,join.fname, folder, data) %>% distinct(),
    features %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))),
    by = c('folder','join.fname')
  ) %>% mutate(iucn_id_no = as.character(iucn_id_no)) %>% dplyr::select(-join.fname,-binomial,-data)
}

# Lastly take the range in case everything fails
if(also_use_eoo && target_range == '_esh'){
  dd <- list.dirs(eoo_path)
  dd <- dd[str_detect(dd,target_resolution)]
  feature_eoo_files <- vector()
  for(folder in dd){
    myLog('Reading in = ', basename(folder))
    ff <- list.files(folder,full.names = TRUE)
    ff <- ff[has_extension(ff,'tif')] 
    feature_eoo_files <- c(feature_eoo_files, ff)
    rm(ff)
  }
  myLog('Loaded extra ', length(feature_eoo_files), ' files')
  features_eoo <- data.frame(fname = basename(feature_eoo_files), full.fname = feature_eoo_files,
                                         folder = basename(dirname(feature_eoo_files)),stringsAsFactors = FALSE) %>% 
    # Get only target resolution
    dplyr::filter(str_detect(folder, target_resolution)) %>% 
    # Subset to only those filenames names in the species_id file
    dplyr::filter(str_to_lower(str_replace_all(fname,' ','_')) %in% str_to_lower(str_replace_all(species_ids$fname,' ','_')) ) %>% 
    dplyr::mutate(type = 'EOO')
  # Join in the id from species id file
  features_eoo <- dplyr::inner_join(
    species_ids %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))) %>% 
      dplyr::select(iucn_id_no, binomial,join.fname, folder, data) %>% distinct(),
    features_eoo %>% mutate(join.fname = str_to_lower(str_replace_all(fname,' ','_'))),
    by = c('folder','join.fname')
  ) %>% mutate(iucn_id_no = as.character(iucn_id_no)) %>% dplyr::select(-join.fname,-binomial,-data)
} else { features_eoo <- data.frame()}

# Create output folder
if(!dir.exists(paste0(fst_path,'/',target_resolution,'_',ifelse(target_range=="",'range','aoh')))){
  dir.create(paste0(fst_path,'/',target_resolution,'_',ifelse(target_range=="",'range','aoh')))
} 
out_dir <- paste0(fst_path,'/',target_resolution,'_',ifelse(target_range=="",'range','aoh')) 

# Ensure that all input feature raster conditions are met
stopifnot(
  assert_that(
    exists('species_ids'),exists('features'),
    exists('features_art_backup'),
    anyDuplicated(features$iucn_id_no) == 0,
    anyNA(features$iucn_id_no) == FALSE,anyNA(features$full.fname) == FALSE,
    res(pu)[1] > 1,
    dir.exists(out_dir),is.number(cores)
  )
)
myLog("All data preparation good. Starting")
# ---------------------------------------------------- #

# Fire up cluster #
cl <- parallel::makeCluster(cores, type = cluster_type, outfile = logfile)
registerDoParallel(cl)

myLog("Preparing data for ", length(pu_cells), " planning units and ", nrow(features) ," remaining features")
# Total length
tl = nrow(features)

# --------------- #
# Run in parallel
o <- foreach(row =  1:tl,
             .inorder = FALSE,
             .errorhandling = 'pass', # SET to stop to check bugs
             .export = c("features","features_art_backup","pu_cells",'pu_id_raster','also_use_eoo','features_eoo',
                         'alignRasters',
                         "species_ids","target_range","out_dir","tl"),
             .packages = c('raster','prioritizr','dplyr','assertthat'),
             .verbose = FALSE) %dopar% {
               # Get the subset
               sub <- features[row,]
               # Define out name
               fname <- paste0(out_dir,'/','id_',unique(sub$iucn_id_no),'__',sub$folder,'__',tools::file_path_sans_ext(sub$fname),'.fst')
               if(file.exists(fname)) {return(NULL)} # skip if output exists already

               myLog("Running ", sub$iucn_id_no ," (",round(row/tl,2)*100,"%)")

               # The feature
               f <- raster::raster( sub$full.fname )
               # Check if projections are matching
               if(compareCRS(f,pu_id_raster)==FALSE){
                 stop('Different projections')
                 print('--> Reprojecting input raster')
                 f <- projectRaster(f,crs = crs(pu_id_raster),method = 'bilinear')
                 f[f==0] <- NA
                 f <- alignRasters(data = f,template = pu_id_raster,method = 'bilinear',func = mean,cl = FALSE)
               }
                 
               # Special case if the AOH is empty of values
               if( (cellStats(f,'max') == 0 || cellStats(f,'max') == -Inf)  && "target_range" != ''){
                  # Use the artifical range instead
                  x <- features_art_backup %>% dplyr::filter(folder == sub$folder, fname == sub$fname)
                  if(nrow(x)==0){ return(NULL)}
                  if(nrow(x)>1){ x <- subset(x, folder == sub$data )}
                  f <- raster::raster( x$full.fname ) # Load the artifical raster instead
                  if(also_use_eoo == FALSE && (cellStats(f,'max') == 0 || cellStats(f,'max') == -Inf) ) { return('Error') }
                  rm(x)
               }
               # If EOO should be used
               if((also_use_eoo && nrow(features_eoo)>0) && (cellStats(f,'max') == 0 || cellStats(f,'max') == -Inf) ){
                 myLog('--> Using EOO instead')
                 sub <- subset(features_eoo, iucn_id_no == sub$iucn_id_no)
                 f <- raster( sub$full.fname )
                 if( (cellStats(f,'max') == 0 || cellStats(f,'max') == -Inf) ) { return('Error') }
               }
               # Check if compare, otherwise align
               if(compareRaster(pu_id_raster,f,stopiffalse = F) == FALSE){
                 f <- alignRasters(data = f,template = pu_id_raster,method = 'bilinear',func = mean,cl = FALSE)
               }
               
               # Security checks
               stopifnot(
                   assert_that(
                    ncell(pu_id_raster) == ncell(f),
                    compareRaster(pu_id_raster,f)
                 )
               )
               # cell indices where feature is present
               cells <- raster::Which(!is.na(f), cells = TRUE)
               # remove cell indices that are not valid planning units, i.e. they have NA costs
               cells <- cells[cells %in% pu_cells]
               # create data.frame with pu ids, feature ids, amounts
               d <- data.frame(pu = pu_id_raster[cells], id = sub$iucn_id_no, amount = f[cells])
               # remove zeros, note we use dplyr because it's actually faster
               d <- dplyr::filter(d, amount > 0)
               if(nrow(d)==0) { return('Error') }
               # Cast file types specifically for saving
               d$pu <- as.integer(d$pu); d$amount <- as.integer(d$amount)
               # Save output
               write.fst(d,  fname)
               return(NULL) # For later
             } # Finish foreach loop
#saveRDS(o,'logs.rds')
myLog('Species without grid cell match: ', length( which(!sapply(o, is.null)) ))
myLog('----------------')
myLog('PU table construction done!')
gc()
#  which(species_ids$id %notin% unique(rij_data$species)) # Manual checking of mismatches
rij <- list.files(out_dir,'*.fst',full.names = T)
myLog('Processed ',length(rij), ' features...' )

# Now construct the global table
# Work around for the larger 10km table
# Split into many smaller chunks
library(data.table);library(progress)
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) # Chunking script

# Make matched data.frame with AOH and range data
if(target_range != ""){
    rij_range <- list.files(paste0(main_path,'results/',target_resolution,'_aoh'),'*.fst',full.names = T)
    ids_range <- str_split( basename(rij_range),'__',simplify = TRUE)[,1]
    # make data.frame of both AOH and range
    rij_df <- data.frame(aoh_id = str_split( basename(rij),'__',simplify = TRUE)[,1], aoh = rij,range = NA,stringsAsFactors = F)
    rij_df$range <- rij_range[match( rij_df$aoh_id, ids_range )]
    rij_df <- rij_df %>% drop_na()
} else {
    rij_df <- data.frame(aoh_id = NA, aoh = NA,range = rij,stringsAsFactors = F)
}

if(run_on_cluster == FALSE | target_resolution == '10km'){
  # Chunck them
  blocks <- chunk(1:nrow(rij_df), 50 )
  # Run in for parallel for both
  rij_data <- data.table()
  pb <- progress_bar$new(total = length(blocks))
  for(ch in 1:length(blocks)){
    x = blocks[[ch]]
    # Process the files in parallel!
    o <- foreach(f =  iter(x),
                 .inorder = FALSE,
                 .errorhandling = 'pass',
                 .export = c('x','rij_df','target_range','pu_data'),
                 .packages = c('fst','dplyr','data.table','assertthat'),
                 .verbose = FALSE) %dopar% {
                   # Load in range
                   range <- read_fst(rij_df$range[f],as.data.table = F)
                   if(target_range != ""){
                     aoh <- read_fst(rij_df$aoh[f],as.data.table = F)
                     assert_that(unique(aoh$id) == unique(range$id))
                     # AOH is by definition nested within range
                     out <- aoh %>% dplyr::filter(pu %in% range$pu)
                   } else {
                     out <- range
                   }
                   # Ensure that pu id's are unique
                   out <- out %>% dplyr::group_by(pu,id) %>% dplyr::summarise(amount = mean(amount)) %>% ungroup()

                   # Convert all pu id's and amounts to integer
                   out$pu <- as.integer(out$pu);out$amount <- as.integer(out$amount)
                   return( out )
                 }
    rij_data <- rbind(rij_data,
                      data.table::rbindlist(o)
    )
    pb$tick()
  }
} else {
  # For 50km data just aggregate them together. Should work?
  rij_data <- data.table::rbindlist(lapply(rij, read.fst))
  # Also ensure that pu_id combinations are unique
  rij_data <- rij_data %>% dplyr::group_by(pu,id) %>% 
    dplyr::summarise(amount = mean(amount)) %>% ungroup()
    
}
rm(o)
# Ensure that all values are between 0 and 1000
rij_data <- subset(rij_data, between(rij_data$amount,0,1000)) 
# - #
myLog("Writing outputs and clean up")

# Write global outputs
write.fst(pu_data,paste0(output_path,'/pu_data_',target_resolution,'.fst') )
writeGeoTiff(pu_id_raster,paste0(output_path,'/puid_',target_resolution,'.tif'),dt =  'INT4U' )

# Sparse matrix
write.fst(rij_data, paste0(output_path,'/',output_filename,'_',target_resolution,target_range,'.fst'), 100)  # use maximum compression and write output
# Clean up #
#lapply(list.files(fst_path,'*.fst',full.names = T), file.remove)
stopImplicitCluster();stopCluster(cl)
gc()
stop('HPC processing DONE!')

#### Load carbon per PU_table ####
# Updated to new carbon data
carbon_soc <- raster(paste0('/media/martin/data/features_esh', "/carbon_agbc/Carbon_SOC_",target_resolution,"_",target_projection,".tif"))
carbon_agbc <- raster(paste0('/media/martin/data/features_esh', "/carbon_agbc/Carbon_AGBC_",target_resolution,"_",target_projection,".tif"))
carbon_bgbc <- raster(paste0('/media/martin/data/features_esh', "/carbon_agbc/Carbon_BGBC_",target_resolution,"_",target_projection,".tif"))
# Align them all
if(!compareRaster(carbon_soc,pu_id_raster,stopiffalse = F)){
  # First for SOC
  carbon_soc <- raster::crop(carbon_soc,  pu_id_raster)
  carbon_soc <- extend(carbon_soc, pu_id_raster)
  carbon_soc <- raster::resample(carbon_soc,pu_id_raster,method = "bilinear")
  # Now for AGBC
  carbon_agbc <- raster::crop(carbon_agbc,  pu_id_raster)
  carbon_agbc <- extend(carbon_agbc, pu_id_raster)
  carbon_agbc <- raster::resample(carbon_agbc,pu_id_raster,method = "bilinear")
  # Now for BGBC
  carbon_bgbc <- raster::crop(carbon_bgbc,  pu_id_raster)
  carbon_bgbc <- extend(carbon_bgbc, pu_id_raster)
  carbon_bgbc <- raster::resample(carbon_bgbc,pu_id_raster,method = "bilinear")
}

# Security checks
stopifnot(
  assert_that(
    compareRaster(carbon_soc,carbon_agbc),
    ncell(pu) == ncell(carbon_soc),
    compareRaster(pu_id_raster,carbon_soc),
    compareRaster(pu_id_raster,carbon_bgbc)
  )
)

## For SOC ##
# cell indices where feature is present
cells_soc <- raster::Which(!is.na(carbon_soc), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells_soc <- cells_soc[cells_soc %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d_soc <- data.frame(pu = pu_id_raster[cells_soc], amount = carbon_soc[cells_soc]) %>% dplyr::filter(amount > 0)
write.fst(d_soc, paste0(output_path,'/pu_carbon_',target_resolution,'_SOC.fst'))

## For AGBC ##
# cell indices where feature is present
cells_agbc <- raster::Which(!is.na(carbon_agbc), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells_agbc <- cells_agbc[cells_agbc %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d_agbc <- data.frame(pu = pu_id_raster[cells_agbc], amount = carbon_agbc[cells_agbc]) %>% dplyr::filter(amount > 0)
write.fst(d_agbc, paste0(output_path,'/pu_carbon_',target_resolution,'_AGBC.fst'))

## For BGBC ##
# cell indices where feature is present
cells_bgbc <- raster::Which(!is.na(carbon_bgbc), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells_bgbc <- cells_bgbc[cells_bgbc %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d_bgbc <- data.frame(pu = pu_id_raster[cells_bgbc], amount = carbon_bgbc[cells_bgbc]) %>% dplyr::filter(amount > 0)
write.fst(d_bgbc, paste0(output_path,'/pu_carbon_',target_resolution,'_BGBC.fst'))

# Combine and aggregate
d <- bind_rows(d_soc,d_agbc,d_bgbc) %>% 
  dplyr::group_by(pu) %>% 
  dplyr::summarise(amount = sum(amount,na.rm = TRUE)) %>% ungroup()

out_name <- paste0(output_path,'/pu_carbon_',target_resolution,'.fst')
write.fst(d, out_name)
rm(d)

#### Load water per PU_table ####
# Updated to new Water data
water <- raster(paste0(biodiversity_path, "/water/cleanwaternorm_pot_",target_resolution,"_",target_projection,".tif"))
# Values are between 0 and 1
# Multiply with 1000 for the prioritization so that amount is comparable to cost
water <- water * 1000
# Align them all
if(!compareRaster(water,pu_id_raster,stopiffalse = F)){
  # First for water
  water <- raster::crop(water,  pu_id_raster)
  water <- extend(water, pu_id_raster)
  water <- raster::resample(water,pu_id_raster,method = "bilinear")
}

# Security checks
stopifnot(
  assert_that(
    ncell(pu) == ncell(water),
    compareRaster(pu_id_raster,water)
  )
)

## For Water ##
# cell indices where feature is present
cells_water <- raster::Which(!is.na(water), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells_water <- cells_water[cells_water %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d_water <- data.frame(pu = pu_id_raster[cells_water], amount = water[cells_water]) %>% dplyr::filter(amount > 0)

out_name <- paste0(output_path,'/pu_water_',target_resolution,'.fst')
write.fst(d_water, out_name)
rm(d_water)

#### Create Protection per PU_table ####
myLog('Extracting Protected Area information per PU')

# Load WDPA data
what <- c('april2019','2010','2014')[1]
wdpa <- raster( paste0(constrain_path, '/wdpa_',what,'_',target_resolution,"_",target_projection,'_',constrain_type,'.tif') )

# Align layers
if(!compareRaster(pu_id_raster, wdpa,stopiffalse = FALSE) | !compareRaster(pu_id_raster, kba,stopiffalse = FALSE)){
  wdpa <- alignRasters(wdpa, pu_id_raster,
                       method = ifelse(constrain_type=="modal","ngb","bilinear"),
                       func = ifelse(constrain_type=="modal",modal,mean) )
}

# Security checks
stopifnot(
  assert_that(
    ncell(pu) == ncell(wdpa),
    compareRaster(pu_id_raster,wdpa)
  )
)

# cell indices where feature is present
cells <- raster::Which(!is.na(wdpa), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells <- cells[cells %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d <- data.frame(pu = pu_id_raster[cells], amount = wdpa[cells]) %>% dplyr::filter(amount > 0)
assert_that(all(between(d$amount,0,1)))
(out_name <- paste0(output_path,'/pu_wdpa_',target_resolution,'_',what,'.fst'))
write.fst(d, out_name)
rm(d)


#### Global biomes and ecoregions ####
# Calculate PU representation per biome and ecoregion
myLog("Processing and assigning Biome and ecoregion data.")
# Load Biome and ecoregion data
biome <- raster( paste0(constrain_path, '/biomes_',target_resolution,
                        ifelse(target_projection=="wgs84","",paste0("_",target_projection)),
                        '_modal.tif') )
ecoregions <- raster( paste0(constrain_path, '/ecoregions_',target_resolution,
                             ifelse(target_projection=="wgs84","",paste0("_",target_projection)),
                             '_modal.tif') )

if(!compareRaster(biome,ecoregions)){
  # Align layers
  biome <- alignRasters(biome, pu_id_raster, method = 'ngb', func = function(x) raster::modal(x, na.rm = TRUE) )
  ecoregions <- alignRasters(ecoregions, pu_id_raster, method = 'ngb', func = function(x) raster::modal(x, na.rm = TRUE) )
}

# Security checks
stopifnot(
  assert_that(
    ncell(pu) == ncell(biome),
    ncell(pu) == ncell(ecoregions),
    compareRaster(pu_id_raster,biome),
    compareRaster(pu_id_raster,ecoregions)
  )
)

# First for biome
# cell indices where feature is present
cells <- raster::Which(!is.na(biome), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells <- cells[cells %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d <- data.frame(pu = pu_id_raster[cells], amount = biome[cells]) %>% dplyr::filter(amount > 0)

out_name <- paste0(output_path,'/pu_biome_',target_resolution,'.fst')
write.fst(d, out_name)
rm(d)

# Then for ecoregion
# cell indices where feature is present
cells <- raster::Which(!is.na(ecoregions), cells = TRUE)
# remove cell indices that are not valid planning units, i.e. they have NA costs
cells <- cells[cells %in% pu_cells]
# create data.frame with pu ids, feature ids, amounts
d <- data.frame(pu = pu_id_raster[cells], amount = ecoregions[cells]) %>% dplyr::filter(amount > 0)

out_name <- paste0(output_path,'/pu_ecoregions_',target_resolution,'.fst')
write.fst(d, out_name)
rm(d)

# ------------------------------------------------- #