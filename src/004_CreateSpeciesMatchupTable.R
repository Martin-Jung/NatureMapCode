library(tidyverse)
library(stringr)
library(lubridate)
library(fst)
library(doParallel)
library(assertthat)
source('src/000_ConvenienceFunctions.R')
cleanNames <- function(x,...){
  # Clean up the biodiversity names 
  g1 <- stringr::str_split(x,"_",simplify = T)[,1]
  g2 <- stringr::str_split(x,"_",simplify = T)[,2]
  g3 <- str_split(stringr::str_split(x,"_",simplify = T)[,3],"\\.",simplify = T)[,1]
  return( data.frame(group = g1, genus = g2, species = g3, fname = x, ..., stringsAsFactors = FALSE) )
}

# Parameters
resolution = "10km" # Target resolution
overwrite = TRUE # Should the table be created anew?

# - Paths - #
output_path       = "data/" # Output path
feature_path      = '/media/martin/data/features/' # Where the features are stored
biodiversity_list = "/media/martin/data/raw/IUCN_RL_2019_2_SpeciesList.csv" # The global biodiversity path lists
gard_list         = "data/gard_ids.rds" # List from Andy Created in Gard processing script
plant_bgci_list   = "data/PLANTS_NON_iucn_gdb.csv" 
plant_miss_list   = "data/plants_wcmc_missing.csv"
if(!file.exists(gard_list)){ stop('Extract shape in GARD EOO script first')}

# ID that is used for creating a plant id (set to one higher than the reptiles)
id_multiplier <- 200000000 # IUCN plant species
# POWO
powo_checklist <- data.table::fread('/media/martin/data/raw/WCSP_KEW/checklist_species.txt') %>% 
  dplyr::filter(acc_full_name_without_family != "Unplaced Unplaced")
powo_checklist_acc <- powo_checklist %>% dplyr::filter(taxon_status_description == 'Accepted')
# Those ones here have duplicated accepted db_id
# powo_checklist_acc[which(powo_checklist_acc$db_id != powo_checklist_acc$accepted_db_id),]

# Biodiversity path(s)
biodiversity_files <- list(
  # GARD data
  reptiles = list.files(paste0(feature_path, "GARDreptilia_",resolution),full.names = F),
  # IUCN reptiles
  iucn_reptiles = list.files(paste0(feature_path, "REPTILIA_",resolution),full.names = F),
  # Birdlife data
  bird = list.files(paste0(feature_path, "AVES_",resolution),full.names = F),
  # Mammal
  mammal = list.files(paste0(feature_path, "MAMMALIA_",resolution),full.names = F),
  # Amphibians
  amphibians = list.files(paste0(feature_path, "AMPHIBIA_",resolution),full.names = F),
  # IUCN plants
  iucn_plantae = list.files(paste0(feature_path, "PLANTAE_",resolution),full.names = F),
  # Kew plants
  kew_iucn_plantae = list.files(paste0(feature_path, "KEW_PLANTAE_IUCN_",resolution),full.names = F),
  kew_bgci_plantae = list.files(paste0(feature_path, "KEW_PLANTAE_BGCI_",resolution),full.names = F),
  # BIEN plants
  BIEN_plantae = list.files(paste0(feature_path, "BIEN_PLANTAE_PPM_",resolution),full.names = F),
  # BIEN plants points
  BIEN_plantae_points = list.files(paste0(feature_path,"BIEN_PLANTAE_POINTS_",resolution),full.names = F),
  # BIEN plants rangebags could not be matched
  # BIEN_plantae_rangebags = list.files(paste0(feature_path, "BIENplantsrangebags_",resolution),full.names = F),
  # New plants species of varying data types
  newPlantPPM = list.files(paste0(feature_path, "newPlantPPM_",resolution),full.names = F),
  newPlantPoints = list.files(paste0(feature_path, "newPlantPoints_",resolution),full.names = F),
  newPlantRangebags = list.files(paste0(feature_path, "newPlantRangebags_",resolution),full.names = F)
) 
# Check that all files have been loaded
assert_that(all( sapply(biodiversity_files, length) >0 ))
myLog("Number species processed: ", sum(sapply(biodiversity_files, length)) )

# ------------------------------------------------------ #
# IUCN data
species <- readr::read_csv(biodiversity_list) %>% 
  dplyr::filter(category %notin% c("EX","EW")) %>% 
  # Now assign
  mutate( fname = paste0(str_to_lower(family),'_',str_replace_all(str_to_lower(binomial)," ","_"),".tif") ) %>% 
  dplyr::select(id_no,binomial,kingdom:criteria,fname) %>% 
  dplyr::rename(iucn_id_no = 'id_no', order = 'order_')
# If Plantae use binomial only to allign with other plant data
species$fname <- ifelse(species$kingdom == 'PLANTAE', paste0( str_replace_all(str_to_lower(species$binomial)," ","_"),".tif"), species$fname)
species_gard <- readr::read_rds(gard_list)
species_missedplants <- read_csv(plant_miss_list)
myLog('Loaded IUCN and GARD data')

# ------------------------------------------- #
#### Process all vertebrates ####
bd_biodiversity1 <- foreach(sn = iter(biodiversity_files$bird),
                       .combine=rbind,
                       .multicombine =TRUE,
                       .inorder = TRUE,
                       .errorhandling = 'pass',
                       .packages = c('tidyverse','stringr'),
                       .verbose = FALSE) %do% {
                         cleanNames(sn, data = 'IUCN Bird')
                       }
bd_biodiversity2 <- foreach(sn = iter(biodiversity_files$reptiles),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'pass',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              cleanNames(sn,data = 'GARD Shai') %>% dplyr::mutate(kingdom = 'REPTILIA')
                            }
bd_biodiversity3 <- foreach(sn = iter(biodiversity_files$mammal),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'pass',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              cleanNames(sn,data = 'IUCN Mammal')
                            }
bd_biodiversity4 <- foreach(sn = iter(biodiversity_files$amphibians),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'pass',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              cleanNames(sn,data = 'IUCN Amphibians')
                            }
bd_biodiversity5 <- foreach(sn = iter(biodiversity_files$iucn_reptiles),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'pass',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              cleanNames(sn,data = 'IUCN Reptiles')
                            }
bd_biodiversity <- bind_rows(bd_biodiversity1,
                             bd_biodiversity2,
                             bd_biodiversity3,
                             bd_biodiversity4,
                             bd_biodiversity5)
rm(bd_biodiversity1,bd_biodiversity2,bd_biodiversity3,bd_biodiversity4,bd_biodiversity5)

myLog('Created matchup table for all included vertebrate species: ', nrow(bd_biodiversity))

# First for all those species that can already be matched
bd_goodmatch <- inner_join(species,bd_biodiversity %>% dplyr::select(fname,data), by = "fname")
# Get the rest
bd_biodiversity <- bd_biodiversity %>% dplyr::filter(fname %notin% bd_goodmatch$fname)

# --- #
# First check those out which are not part of the GARD dataset
bd_others  <- bd_biodiversity %>% dplyr::filter(data != "GARD Shai") %>% 
  # Correct genus heading and remove irrelevant columns
  mutate(genus = str_to_title(genus)) %>% dplyr::select(-group)

if(nrow(bd_others) > 0){
  # Check a specific known mismatches with species
  x1 <- species %>% dplyr::filter(binomial %in% c('Alouatta seniculus ssp. juara','Oryx beisa ssp. beisa','Lyciasalamandra billae ssp. billae') ) %>%
    # Get species name
    mutate(binomial = paste0(str_split(binomial," ",simplify = T)[,1]," ",str_split(binomial," ",simplify = T)[,2]) )
  x1$data <- c('IUCN Mammal','IUCN Mammal','IUCN Amphibians')
  if(nrow(bd_others)>0){
    bd_others <- left_join(bd_others,x1, by = c("genus","species"))
  } else {
    bd_others <- x1
  }
  rm(x1)
  # Apparend to goodmatch
  bd_goodmatch <- bind_rows(bd_goodmatch,bd_others)
}
rm(bd_others)

# Finally do the reptiles from GARD et al.
bd_reptiles <- bd_biodiversity %>% dplyr::filter(data == "GARD Shai") %>% 
  # Make title
  mutate(genus = str_to_title(genus)) %>% 
  dplyr::select(-group) %>% 
  # Construct missing binomial 
  mutate( binomial = paste0(genus," ", species))

assert_that(all(bd_reptiles$binomial %in% species_gard$binomial) )
# Merge with species gard information
bd_reptiles <- left_join(
  bd_reptiles %>% dplyr::select(genus, species, fname, data, binomial),
  species_gard %>% dplyr::select(gard_id,class,family,order,binomial) %>% 
    dplyr::mutate(kingdom = 'ANIMALIA') %>% 
    # Use gard_id in case there is a missing iucn_id_no
    dplyr::rename(iucn_id_no = gard_id),
  by = 'binomial'
)
myLog("Now still ", length(which(is.na(bd_reptiles$iucn_id_no))), " reptiles without IUCN / Gard ID" )

bd_goodmatch <- bind_rows(bd_goodmatch,bd_reptiles)
rm(bd_reptiles,bd_biodiversity) # Clean

# Remove eventual duplicates
bd_goodmatch <- bd_goodmatch %>% distinct() %>% 
  # Remove species column (creates problems later with prioritization). Can be recreated from filename 
  dplyr::select(-species)

# Check again that everything is there
stopifnot(
  assert_that(anyNA(bd_goodmatch$iucn_id_no) == FALSE,
              anyNA(bd_goodmatch$binomial) == FALSE,
              anyNA(bd_goodmatch$kingdom) == FALSE)
)
# Correct some misclassifications (from IUCN script)
# FIXME: This should ideally be done in the GARD script 
bd_goodmatch[which(bd_goodmatch$class=='Insecta'),c('class','order','family')] <- c('Reptilia','Squamata','Dipsadidae')
bd_goodmatch$class[bd_goodmatch$class=='Reptilia'] <- "REPTILIA"
bd_goodmatch$phylum[is.na(bd_goodmatch$phylum)] <- "CHORDATA"
bd_goodmatch$order[!is.na(bd_goodmatch$order)] <- str_to_upper(bd_goodmatch$order[!is.na(bd_goodmatch$order)])
bd_goodmatch$family[!is.na(bd_goodmatch$family)] <- str_to_upper(bd_goodmatch$family[!is.na(bd_goodmatch$family)])

# Make column to link to the folder directory
bd_goodmatch$folder <- paste0( str_to_upper(bd_goodmatch$class),"_",resolution )
bd_goodmatch$folder[which(bd_goodmatch$data=="GARD Shai")] <- paste0('GARDreptilia_',resolution) # Different folder for IUCN reptiles

# Final Security checks
assertthat::assert_that(
  any(!is.na(bd_goodmatch$iucn_id_no)),
  any(!duplicated(bd_goodmatch$fname)),
  all(bd_goodmatch$fname %in% basename(unlist(biodiversity_files) ))
)

#### Adding plants ####
# Process the plant names and assign an id to each of it.
# Start formatting and processing
bd_biodiversity1 <- foreach(sn = iter(biodiversity_files$iucn_plantae),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                data.frame(
                                  group = NA,
                                  kingdom = 'PLANTAE',
                                  genus  = stringr::str_split(sn,"_",simplify = T)[,1],
                                  species = tools::file_path_sans_ext(stringr::str_split(sn,"_",simplify = T)[,2]),
                                  fname = sn,
                                  data = 'IUCN Plants',
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
# Add Folder
bd_biodiversity1$folder <- paste0("PLANTAE_",resolution)
assert_that(
  anyNA(bd_biodiversity1$genus) == FALSE,
  anyNA(bd_biodiversity1$species) == FALSE
)
# --- #
# Now all the rest
# KEW - IUCN
bd_biodiversity2 <- foreach(sn = iter(biodiversity_files$kew_iucn_plantae),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                data.frame(
                                  group = NA,
                                  kingdom = 'PLANTAE',
                                  genus  = stringr::str_split(sn," ",simplify = T)[,1],
                                  species = tools::file_path_sans_ext(stringr::str_split(sn," ",simplify = T)[,2]),
                                  fname = sn,
                                  data = 'Kew IUCN Plants',
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
# Add Folder
bd_biodiversity2$folder <- paste0("KEW_PLANTAE_IUCN_",resolution)
assert_that(
  anyNA(bd_biodiversity2$genus) == FALSE,
  anyNA(bd_biodiversity2$species) == FALSE
)

# KEW - BGCI
biodiversity_files$kew_bgci_plantae <- str_replace(biodiversity_files$kew_bgci_plantae,'_',' ') # 5 species or so with underscore!
bd_biodiversity3 <- foreach(sn = iter(biodiversity_files$kew_bgci_plantae),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                data.frame(
                                  group = NA,
                                  kingdom = 'PLANTAE',
                                  genus  = stringr::str_split(sn," ",simplify = T)[,1],
                                  species = tools::file_path_sans_ext(stringr::str_split(sn," ",simplify = T)[,2]),
                                  fname = sn,
                                  data = 'Kew BGCI Plants',
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
# Add Folder
bd_biodiversity3$folder <- paste0("KEW_PLANTAE_BGCI_",resolution)
assert_that(
  anyNA(bd_biodiversity3$genus) == FALSE,
  anyNA(bd_biodiversity3$species) == FALSE
)

bd_biodiversity4 <- foreach(sn = iter(biodiversity_files$BIEN_plantae),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                    data.frame(
                                      group = NA,
                                      kingdom = 'PLANTAE',
                                      genus  = stringr::str_split(str_replace_all(sn," ","_"),"_",simplify = T)[,1],
                                      species = tools::file_path_sans_ext(stringr::str_split(str_replace_all(sn," ","_"),"_",simplify = T)[,2]),
                                      fname = sn,
                                      data = 'BIEN Plantae PPM',
                                      stringsAsFactors = FALSE
                                    )
                              )
                            }
# Add Folder
bd_biodiversity4$folder <- paste0("BIEN_PLANTAE_PPM_",resolution)
assert_that(
  anyNA(bd_biodiversity4$genus) == FALSE,
  anyNA(bd_biodiversity4$species) == FALSE
)

bd_biodiversity5 <- foreach(sn = iter(biodiversity_files$BIEN_plantae_points),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                data.frame(
                                  group = NA,
                                  kingdom = 'PLANTAE',
                                  genus  = stringr::str_split(str_replace_all(sn," ","_"),"_",simplify = T)[,1],
                                  species = tools::file_path_sans_ext(stringr::str_split(str_replace_all(sn," ","_"),"_",simplify = T)[,2]),
                                  fname = sn,
                                  data = 'BIEN Plantae Points',
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
# Add Folder
bd_biodiversity5$folder <- paste0("BIEN_PLANTAE_POINTS_",resolution)
assert_that(
  anyNA(bd_biodiversity5$genus) == FALSE,
  anyNA(bd_biodiversity5$species) == FALSE
)

# -- Our new plant species addded -- #
bd_biodiversity6 <- foreach(sn = iter(biodiversity_files$newPlantPPM),
                            .combine=rbind,
                            .multicombine =TRUE,
                            .inorder = TRUE,
                            .errorhandling = 'stop',
                            .packages = c('tidyverse','stringr'),
                            .verbose = FALSE) %do% {
                              return(
                                data.frame(
                                  group = NA,
                                  kingdom = 'PLANTAE',
                                  genus  = stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,1],
                                  species = tools::file_path_sans_ext(stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,2]),
                                  fname = sn,
                                  data = 'New Plants - PPM',
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
# Add Folder
bd_biodiversity6$folder <- paste0("newPlantPPM_",resolution)
assert_that(
  anyNA(bd_biodiversity6$genus) == FALSE,
  anyNA(bd_biodiversity6$species) == FALSE
)
bd_biodiversity7 <- foreach(sn = iter(biodiversity_files$newPlantPoints),
                             .combine=rbind,
                             .multicombine =TRUE,
                             .inorder = TRUE,
                             .errorhandling = 'stop',
                             .packages = c('tidyverse','stringr'),
                             .verbose = FALSE) %do% {
                               return(
                                 data.frame(
                                   group = NA,
                                   kingdom = 'PLANTAE',
                                   genus  = stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,1],
                                   species = tools::file_path_sans_ext(stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,2]),
                                   fname = sn,
                                   data = 'New Plants - Points',
                                   stringsAsFactors = FALSE
                                 )
                               )
                             }
# Add Folder
bd_biodiversity7$folder <- paste0("newPlantPoints_",resolution)
assert_that(
  anyNA(bd_biodiversity7$genus) == FALSE,
  anyNA(bd_biodiversity7$species) == FALSE
)
bd_biodiversity8 <- foreach(sn = iter(biodiversity_files$newPlantRangebags),
                             .combine=rbind,
                             .multicombine =TRUE,
                             .inorder = TRUE,
                             .errorhandling = 'stop',
                             .packages = c('tidyverse','stringr'),
                             .verbose = FALSE) %do% {
                               return(
                                 data.frame(
                                   group = NA,
                                   genus  = stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,1],
                                   species = tools::file_path_sans_ext(stringr::str_split(str_replace(sn, ' ','_'),"_",simplify = T)[,2]),
                                   fname = sn,
                                   data = 'New Plants - Rangebag',
                                   stringsAsFactors = FALSE
                                 )
                               )
                             }
# Add Folder
bd_biodiversity8$folder <- paste0("newPlantRangebags_",resolution)
assert_that(
  anyNA(bd_biodiversity8$genus) == FALSE,
  anyNA(bd_biodiversity8$species) == FALSE
)
# ----------------------------------------- #
# For the plants in order:
# Process first the IUCN / KEW and BGCI as expert data
# Then the new models of species
# Then everything else not already in the new species
# -> All duplicates are inserted Do the actual data selection in the PU table script

# First for all those species that can already be matched
plants_goodmatch <- inner_join(species, bd_biodiversity1 %>% dplyr::select(fname,data,folder), by = "fname")

# New species from KEW
plants_kew <- bind_rows(
  bd_biodiversity2,
  bd_biodiversity3
) %>% dplyr::mutate(genus = str_to_sentence(genus), species = str_to_lower(species)) %>% # Correct matching names
  dplyr::mutate( binomial  = paste0(str_to_sentence(genus)," ", str_to_lower(species) ) ) %>% 
  dplyr::distinct() # Unique values

# Load in the plant IDs from WCMC and join them back in
new_plant_ids <- read_csv( plant_bgci_list ) %>% dplyr::mutate(BINOMIAL = str_to_sentence( str_trim(BINOMIAL,side = 'both') ) )
names(new_plant_ids) <- str_to_lower(names(new_plant_ids)) # lower case column names
# Join in
plants_kew <- left_join(
  plants_kew %>% dplyr::select(binomial,data,folder,fname),
  new_plant_ids,
  by = 'binomial'
) %>% dplyr::rename(iucn_id_no = id_no)
rm(new_plant_ids)

# Check for all species that are already in the good matching ones and insert the missing
plants_goodmatch <- bind_rows(
  plants_goodmatch,
  plants_kew#[which(plants_kew$binomial %notin% str_to_sentence(plants_goodmatch$binomial)),]
)
rm(plants_kew)
assert_that( anyNA(plants_goodmatch$binomial)==FALSE, anyNA(plants_goodmatch$fname) == FALSE )

# Then the newly modelled global species
plants_global <- bind_rows(
  bd_biodiversity6,
  bd_biodiversity7,
  bd_biodiversity8
) %>% dplyr::mutate(genus = str_to_sentence(genus), species = str_to_lower(species)) %>% # Correct matching names
  dplyr::mutate( binomial  = paste0(str_to_sentence(genus)," ", str_to_lower(species) ) ) %>% 
  dplyr::distinct() # Unique values

# Check for all species that are already in the good matching ones and insert the missing
plants_goodmatch <- bind_rows(
  plants_goodmatch,
  plants_global#[which(plants_global$binomial %notin% str_to_lower(plants_goodmatch$binomial)),]
)
rm(plants_global)
assert_that( anyNA(plants_goodmatch$binomial)==FALSE, anyNA(plants_goodmatch$fname) == FALSE )

# Then get the rest from BIEN
plants_bien <- bind_rows(bd_biodiversity4,
                         bd_biodiversity5
                             ) %>% 
  dplyr::mutate(genus = str_to_sentence(genus), species = str_to_lower(species)) %>% # Correct matching names
  dplyr::mutate( binomial  = paste0(str_to_sentence(genus)," ", str_to_lower(species) ) ) %>% 
  dplyr::distinct() # Unique values

# Check for all species that are already in the good matching ones and insert the missing
plants_goodmatch <- bind_rows(
  plants_goodmatch,
  plants_bien#[which(plants_bien$binomial %notin% str_to_lower(plants_goodmatch$binomial)),]
)
rm(plants_bien)
assert_that( anyNA(plants_goodmatch$binomial)==FALSE, anyNA(plants_goodmatch$fname) == FALSE )

rm(bd_biodiversity1,bd_biodiversity2,bd_biodiversity3,bd_biodiversity4,bd_biodiversity5,bd_biodiversity6,bd_biodiversity7,bd_biodiversity8)

# Remove files without fname
plants_goodmatch <- plants_goodmatch %>% drop_na(fname)

# Last additions
plants_goodmatch$kingdom <- ifelse(is.na(plants_goodmatch$kingdom),'PLANTAE',plants_goodmatch$kingdom)
plants_goodmatch$group <- NULL
plants_goodmatch$species <- NULL

# Some checks again
assert_that(
  anyNA(plants_goodmatch$binomial)==FALSE,
  anyNA(plants_goodmatch$folder)==FALSE,
  anyNA(plants_goodmatch$fname)==FALSE,
  all(!duplicated(paste0(plants_goodmatch$folder,'/',plants_goodmatch$fname)) )
)
# ------------------------------ #
# Now merge in the POWO ID to the accepted species names
plant_ids <- powo_checklist_acc %>% dplyr::rename(binomial = acc_full_name_without_family) %>% dplyr::select(binomial, db_id) %>%
  dplyr::mutate(binomial = str_to_sentence(binomial)) %>% distinct()
# For those duplicated binomial names, take the first name/ID combination
plant_ids <- plant_ids[-which(duplicated(plant_ids$binomial)),]
assert_that( anyDuplicated(plant_ids$binomial) == 0, anyDuplicated(plant_ids$db_id) == 0 )

# Join in
plants_goodmatch <- dplyr::left_join(
 plants_goodmatch %>% dplyr::mutate(binomial = str_to_sentence(binomial) ),
 plant_ids,
 by = 'binomial'
      );rm(plant_ids)
# Check also for the synonym's in POWO
plants_noid <- plants_goodmatch %>% dplyr::filter(is.na(db_id)) %>% dplyr::select(-db_id)
plant_ids <- powo_checklist %>% dplyr::select(full_name_without_family,acc_full_name_without_family, db_id) %>% distinct()
plant_ids <- plant_ids[-which(duplicated(plant_ids$full_name_without_family)),]
assert_that( anyDuplicated(plant_ids$full_name_without_family) == 0, anyDuplicated(plant_ids$db_id) == 0 )
plants_noid <- left_join(plants_noid, plant_ids, by = c('binomial' = 'full_name_without_family'))
# --- #
plants_goodmatch <- bind_rows(
  plants_goodmatch  %>% dplyr::filter(!is.na(db_id)),
  plants_noid %>% dplyr::select(-acc_full_name_without_family)
)
myLog('What proportion not matched to POWO? : ', round(length(which(is.na(plants_goodmatch$db_id))) / nrow(plants_goodmatch),3) )
rm(plant_ids,plants_noid)
# ---------------------------- #
# Assign new ids
plants_withid <- plants_goodmatch %>% dplyr::filter(!is.na(iucn_id_no)) %>% distinct()
plants_without <- plants_goodmatch %>% dplyr::filter(is.na(iucn_id_no)) %>% dplyr::select(-iucn_id_no) %>% distinct()
# Assign same species id to those matching already in there
pid <- plants_withid %>% dplyr::select(binomial,iucn_id_no) %>% distinct()
pid <- pid[-which(duplicated(pid$binomial)),]
plants_without <- left_join(plants_without, pid, by = 'binomial');rm(pid)
plants_withid <- bind_rows(plants_withid, plants_without %>% dplyr::filter(!is.na(iucn_id_no)) )
plants_without <- plants_without %>% dplyr::filter(is.na(iucn_id_no)) %>% dplyr::select(-iucn_id_no)
assert_that((nrow(plants_withid) + nrow(plants_without)) == nrow(plants_goodmatch))
# Now for the rest
# Make a unique new id based on the species
plants_without <- left_join(plants_without,
                            plants_without %>% dplyr::select(binomial) %>% distinct() %>% 
                              dplyr::mutate(iucn_id_no = row_number() + id_multiplier ),
                            by = 'binomial')
# Security check
assert_that(anyNA(plants_without$iucn_id_no) == FALSE, anyNA(plants_withid$iucn_id_no) == FALSE,
            (nrow(plants_withid) + nrow(plants_without)) == nrow(plants_goodmatch) )
plants_goodmatch <- bind_rows(plants_withid, plants_without) # now merge them back together
rm(plants_withid,plants_without)

# Security checks
assert_that(
  any(!is.na(plants_goodmatch$iucn_id_no)),
  any(!is.na(plants_goodmatch$folder)),
  any(!duplicated(plants_goodmatch$fname)),
  any(duplicated(paste0(plants_goodmatch$folder,'/',plants_goodmatch$fname))) == FALSE,
  all( names(bd_goodmatch) %in% names(plants_goodmatch) )
)
myLog('Created matchup table for all included species: ', nrow(plants_goodmatch))

bd_goodmatch <- bind_rows(
  bd_goodmatch,
  plants_goodmatch
)

# ----------------------------------- #
# Save the final output 
# ----------------------------------- #
write.fst(bd_goodmatch,paste0(output_path,"speciesID_table_",resolution,".fst"))
rm(bd_goodmatch)
