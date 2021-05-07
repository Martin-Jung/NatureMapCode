#### Start of parameter setting ####
# This script specifies and solves all global conservation problems
# A number of different problems are formulated, namely
if(dir.exists('/h/u046/jung/naturemap')) { setwd('/h/u046/jung/naturemap');run_on_cluster <- TRUE } else { run_on_cluster <- FALSE }
# Number of cores and cluster_type for parallel processing 
cores <- ifelse(run_on_cluster, 1, 7 )
cores_repr <- ifelse(run_on_cluster, 1, 7 ) # The number of cores for the representation calculation
set.seed(3663) # set seed for reproducibility
overwrite_existing_data <- FALSE # If set to true, don't skip existing data

# Paths
output_path <- "data" # Output path for database  
target_resolution <- "10km" # Target resolution
target_resolution_number <- as.numeric(gsub("\\D","",target_resolution)) # Target resolution as number 
target_range <- c("","_esh")[2]

# RIJ table   
rij_path <- paste0(output_path,"/table_puspecies",'_',target_resolution,target_range,'.fst') # The final output table

# PU data
pu_path <- paste0(output_path,'/pu_data','_',target_resolution,'.fst')
puid_path <- paste0(output_path,'/puid','_',target_resolution,'.tif')

# Carbon data
carbon_path <- paste0(output_path,"/pu_carbon_",target_resolution,".fst")
# Water path
water_path <- paste0(output_path,"/pu_water_",target_resolution,".fst")

# Feature data
feature_path <- paste0(output_path,"/speciesID_table_",target_resolution,".fst")
constraints_path <- c(
  paste0(output_path,"/pu_wdpa_",target_resolution,"_april2019.fst")
)
biome_path <- paste0(output_path,"/pu_biome_",target_resolution,".fst")
ecoregion_path <- paste0(output_path,"/pu_ecoregions_",target_resolution,".fst")
phylo_path <- paste0(output_path,"/ED_threat_scores.fst")

# Results path
results_path <- paste0('results/',target_resolution,target_range)
if(!dir.exists(results_path)){ dir.create(results_path) }

# ---- Prioritize targets ---- #
# Maxfeat_speciestargets calculates the shortfall. The others should not be used.
what <- c('minshort_speciestargets','maxutilityobjective','maxcoverobjective','datapreparation')[1]
target_datasets <- c("IUCN Mammal","IUCN Reptiles","IUCN Bird","IUCN Amphibians","GARD Shai") # Vertebrate datasets
# Biodiversity parameters
exclude_plants <- FALSE # Don't include plants
phylo_weights <- FALSE # Phylogenetic weights by species
phylo_comparison <- FALSE # Make a comparison by subsetting only to those species with a score
phylo_type <- c('rank_weight','ed','edge','noweight')[1] # Which species weight specifically ?
# Other data sources  
include_carbon <- FALSE # Should carbon be included as a feature? 
include_water <- FALSE # Should carbon be inclu ded as a feature? 
carbon_multiplier <- 1 # Should the carbon weight be reduced?
water_multiplier <- 1 # Should the water weight be reduced?       
  
# -- ID subsetting --             
repr_id <- TRUE # Should there be a subsetting  
repr_id_path <- 'sets' # Path to representative IDs
repr_numbers <- 1:10 #1  Which representative sets to run  
# -------- #
  
## Parameters for prioritization ##     
budgets <- seq(0.1,1,0.1) # Set budgets 
decision_type <- c("binary","proportion")[2] # What kind of decision type?
# Varieties
split_species <- TRUE # Should species ranges be split ?
split_id <- c("biome.id","ecoregion.id")[1] # Which id used to split ?  
pa_lockedin <- TRUE # Should protected areas be locked in?    
ensure_all_land <- FALSE # Include the pu_data as feature to ensure that all features are selected
port <- FALSE # Calculate portfolios
run_variety <- c('','carbonweights','carbonwaterweights','carbononly')[1] # Should a different run be conducted?
calculate_representation <- TRUE # Calculate the representation
optimality_gap <- 0 # The gap at which the solver can stop (0 = Optimality)
time_limit <- .Machine$integer.max # Should the computation be terminated at some point?

# -----------------------------  #
# End of parameter setting here  #
# -----------------------------  #    
# load packages
library(prioritizr)
library(data.table)
library(raster)
library(parallel)
library(doParallel)
library(scales)
library(fst)
library(assertthat)
library(tidyverse);require(progress)
source('src/000_ConvenienceFunctions.R')
# ------------------- #
# Security checks
stopifnot(
  assert_that(
    file.exists(pu_path),
    file.exists(carbon_path),
    file.exists(water_path),
    file.exists(feature_path),
    all(sapply(constraints_path, file.exists)),
    file.exists(biome_path),
    file.exists(phylo_path),
    file.exists(ecoregion_path),
    is.number(carbon_multiplier),is.number(water_multiplier)
  )
)
if(pa_lockedin & decision_type == 'binary'){ warning('Protected area lock in and binary decisions unlikely to work...') }
if(split_id == "ecoregion.id") { stop('Ecoregion id split Not yet implemented!!!')}
if(phylo_weights & (exclude_plants == FALSE & phylo_type != 'rank_weight')) { warning('Phylogenetic weights for plants not incorporated. Removing plants');exclude_plants <- TRUE }
if(phylo_weights & !repr_id) { warning('Run phylogenetic weights on all vertebrates');repr_id <- FALSE }
if(run_variety == 'carbononly' & include_carbon == FALSE) { stop('This variety does not make sense unless carbon is included') }
if(repr_id && repr_id_path == '') { stop('Set a path for the representative ids') }
# ---------------- #
if(repr_id) { 
  lf <- list.files(repr_id_path,'splist',full.names = T)
  repr_id_file <- lapply(lf, function(x) { read_csv(x) %>% dplyr::select(2) } )
  names(repr_id_file) <- paste0('id',1:length(repr_id_file))
  # Rather quick and dirty rename the columns too
  for(i in names(repr_id_file)){
    x = repr_id_file[[i]]; names(x) <- i
    repr_id_file[[i]] <- x
  }
  rm(lf)
  if(run_variety!=''){
    # If a run_variety is species, use a combined three sets
    repr_id_file2 <- list()
    repr_id_file2[['id1']] <- data.frame(id1= unique(c( repr_id_file[[1]][['id1']] , repr_id_file[[2]][['id2']], repr_id_file[[3]][['id3']])) )
    repr_id_file <- repr_id_file2;rm(repr_id_file2)
  }
} else {repr_id_file <- data.frame(normalrun = NA);repr_numbers <- 1 }
#permut_run = names(repr_id_file)[1]  
for(permut_run in names(repr_id_file)[repr_numbers] ){
myLog('Doing a representative run for species set ', permut_run)

# Load in data
myLog("...Reading in data...", target_range, " at ", target_resolution)
pu_data <- read_fst(pu_path)
carbon_data <- read_fst(carbon_path) %>%
  mutate(id = 9999999, species = 9999999) %>% # Assign ID to carbon data
  mutate(amount = amount * 100) %>% # Convert ha to squarekm2
  # Left_join in carbon with fractional cost
  dplyr::left_join(., pu_data %>% dplyr::mutate(cost = cost / 1000), by = c('pu' = 'id') ) %>% 
  mutate(amount = amount * cost) %>% dplyr::select(-cost)

if(include_water){
  water_data <- read_fst(water_path) %>% 
    mutate(id = 9999998, species = 9999998) %>% 
    # Also multiply water with fractional cost
    dplyr::left_join(., pu_data %>% dplyr::mutate(cost = cost / 1000), by = c('pu' = 'id') ) %>% 
    mutate(amount = amount * cost) %>% dplyr::select(-cost)
}
# Load PU id raster
pu_id_raster <- raster(puid_path)
pu_id_raster[pu_id_raster==0] <- NA

# Regions data
regions_data <- full_join(
  read_fst(biome_path) %>% dplyr::rename(biome.id = 'amount',id = 'pu'),
  read_fst(ecoregion_path) %>% dplyr::rename(ecoregion.id = 'amount',id = 'pu'),
  by = 'id'
) %>% drop_na() # Remove missing data

# -- Prepare feature input data -- ## 
# Load feature data
feature_data <- read_fst(feature_path) %>% dplyr::rename(id = iucn_id_no)
stopifnot( assert_that(!anyNA(feature_data$binomial)) )
# Always take IUCN before GARD
ids <- feature_data$id[which(feature_data$data == "GARD Shai")][feature_data$binomial[which(feature_data$data == "GARD Shai")] %in% feature_data$binomial[which(feature_data$data == "IUCN Reptiles")] ]
feature_data <- feature_data %>% dplyr::filter(id %notin% ids)
rm(ids)

# No duplicate of species
if(exclude_plants==FALSE){
  plant_ids <- feature_data %>% dplyr::filter(kingdom == 'PLANTAE')
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
  plant_ids3 <- bind_rows(plant_ids3_1,plant_ids3_2)
  
  plant_ids <- bind_rows(plant_ids1,plant_ids2,plant_ids3) %>% distinct();rm(plant_ids1,plant_ids2,plant_ids3)
  stopifnot( assert_that(anyDuplicated(plant_ids$id) == 0,
                         anyNA(plant_ids$id) == FALSE )
  )
  # Get the IDs for the feature dataset
  ids <- c( feature_data$id[which(feature_data$id %in%  plant_ids$id)],
            feature_data$id[which(feature_data$data %in% target_datasets)]
  )
  rm(plant_ids)
} else{
  ids <- c(feature_data$id[which(feature_data$data %in% target_datasets)] )
}
stopifnot( assert_that(!anyNA(ids)) ) # Sec check 

# ---------------------------#
# From here process only representative sets
if(repr_id){
  # Load the subsampled estimates
  subsample <- repr_id_file[[permut_run]] # Load the subsample 
  ids  <- ids[ids %in% subsample[[permut_run]] ] # Remove all ids not part of the representative set
  stopifnot(assert_that(length(ids) > 0 ))
}

# Filter feature data
feature_data <- feature_data %>% dplyr::filter(id %in% ids)

# -------- Rij table load in  ------- #
stopifnot( assert_that(file.exists(rij_path)) )
# Load in PU table directly
rij_data <- read_fst(rij_path,as.data.table = TRUE)
rij_data <- subset(rij_data, id %in% ids) # First subset
# ----RIJ table end---- #
rm(ids)

# The species AOH amount can never exceed the budget 
# Ensure that this condition is met
rij_data <- rij_data %>% 
  left_join(., pu_data, by = c("pu" = "id") ) %>% 
  dplyr::mutate( amount = pmin(amount,cost) ) %>% 
  dplyr::select( - cost )

# Security checks
stopifnot(
  assert_that(
    res(pu_id_raster)[1]*1e-3 == target_resolution_number,  # Check that resolution is identical to pu cell size resolution in meters
    cellStats(pu_id_raster,"max") == n_distinct(pu_data$id),# Number of pu units correct?
    max(rij_data$amount) <= max(pu_data$cost), # The amount cannot be greater than the cost!
    all(unique(rij_data$id) %in% feature_data$id) # All species present in feature data?
  )
)
suppressMessages(gc(verbose = FALSE) )

#### ---- Target setting ---- ####
myLog("Preparing targets")
# Should species
if(split_species){
  myLog("Splitting by target group: ", split_id, "- Number distinct species in rij_table: ", length(unique(rij_data$id)) )
  
  # Append region data with biomes and ecoregions
  rij_data <- dplyr::left_join(rij_data %>% dplyr::mutate(pu = as.character(pu)),
                               regions_data %>% dplyr::mutate(id = as.character(id)), by = c("pu" = "id")) %>% 
    dplyr::mutate(
      fid = id, species = paste0(id,"__",biome.id)
    )
  rij_data <- rij_data %>% drop_na(all_of(split_id))

  # Now calculate targets
  rij_data$area <- (target_resolution_number * target_resolution_number) # Save an area amount estimate
  
  # --- #
  # Criteria being that at least 2200km are part of the new biome area. If not. assign area to the closest, contigious biome
  problematic_species <- rij_data %>% dplyr::group_by(id,species) %>% 
    dplyr::summarise(feature_abundance = sum((amount / 1000) * area) ) %>% 
    # Find those problematic new species splits
    dplyr::filter(feature_abundance < 2200) %>% dplyr::select(id,species) %>% ungroup() %>% distinct() %>% 
    dplyr::rename(problematic_species = species)
  # And find the largest id's
  largest_species <- rij_data %>% dplyr::select(id,species,amount,area) %>% 
    dplyr::group_by(id,species) %>% 
    dplyr::summarise(feature_abundance = sum((amount / 1000) * area) ) %>% 
    dplyr::top_n(1, feature_abundance ) %>% ungroup() %>% dplyr::select(id,species) %>% distinct() 
  # Remove id's in problematic species and join them to the largest one
  # E.g. left_join prob
  species_match <- left_join(largest_species, problematic_species, by = 'id') %>% 
    drop_na() %>% # Drop non-problematic species
    dplyr::rename(new_species = species) %>% dplyr::select(-id) 
  # All species in problematic_species should be set to species
  rij_data <- rij_data %>%
    dplyr::left_join(., species_match, by = c('species' = 'problematic_species')) %>% 
    # If non
    dplyr::mutate(species = ifelse(is.na(new_species), species, new_species) ) %>% 
    dplyr::select(-new_species)
  rm(species_match, problematic_species, largest_species ) # Cleanup
  # --- #
  # Calculate targets
  species_targets <- calc_targets2(rij_data,multiplier = 1000) %>% mutate(species = factor(species))
  # Join in original species_id
  species_targets <- left_join(species_targets,rij_data %>% dplyr::select(species,fid) %>% distinct(),
                               by = "species") 
  # Convert fid to numeric
  species_targets$fid <- as.numeric(as.character(species_targets$fid))
  
  stopifnot(
    assert_that(
      any(duplicated(species_targets$species)) == FALSE,
      length(which(species_targets$relative_target>1)) == 0,
      n_distinct(rij_data$fid) == n_distinct(species_targets$fid),
      n_distinct(species_targets$species) == n_distinct(rij_data$species),
      anyNA(species_targets$fid) == FALSE,
      all( species_targets$fid == str_split(species_targets$species,"__",simplify = T)[,1]  )
    )
  )
  # Need to rename otherwise throws an error
  feature_data <- feature_data %>%
    # Use the feature id for joining
    dplyr::mutate(fid = as.numeric( as.character(id) ) ) %>% dplyr::select(-id) %>% 
    dplyr::mutate(fname = str_replace(fname," ","_") ) %>% 
    dplyr::select(fid,binomial) %>% distinct() %>% # 
    left_join(., species_targets,by = "fid") %>% 
    # Fill missing targets if species not represented in pu table 
    mutate(relative_target = ifelse(is.na(relative_target),0,relative_target) )
  # Remove species not represented in PU table or that have no target
  feature_data <- feature_data[-which(is.na(feature_data$species)),]
  stopifnot( assert_that( length( which( duplicated( feature_data$species ))) == 0) ) # Assess this
  
  # Set and correct names to be the file name
  feature_data$name <- feature_data$binomial
  feature_data$name <- paste0(feature_data$name,"__",str_split(feature_data$species,"__",simplify = T)[,2]) # Correct file name
  
  # Number at high numeric ranges are not stored correctly in R. Therefore reset the counter to 1:nrow()
  new_ids <- data.frame(species = feature_data$species, newid = 1:nrow(feature_data))
  # Check whether these artifical ids are not already existing as feature_id
  stopifnot(
    assert_that(
      anyDuplicated(new_ids$species) ==FALSE,
      any(new_ids %in% as.numeric( as.character(feature_data$fid) ) ) == FALSE
    )
  )
  feature_data$id <- new_ids$newid
  rij_data <- rij_data %>% dplyr::left_join(., new_ids, by = 'species') %>% 
    dplyr::select(-species) %>% dplyr::rename(species = "newid")
  
  # Further assertions
  stopifnot( assert_that( all( feature_data$id %in% unique(rij_data$species) ) ) )
  
  # Clean up and keep only those columns necessary for the prioritization
  rij_data <- rij_data %>% dplyr::select(pu,species,amount) %>% #dplyr::rename(id = "fid") %>% 
    dplyr::mutate(pu = as.numeric(as.character(pu)), species = as.numeric(as.character(species))) #,id = as.numeric(id) )
  feature_data <- feature_data %>% dplyr::select(fid,id,name,relative_target,species) %>% distinct()
} else 
  {
  rij_data$species <- paste0(rij_data$id)
  rij_data$area <- (target_resolution_number * target_resolution_number) # Save an area amount estimate
  species_targets <- calc_targets2(rij_data,multiplier = 1000) %>% mutate(species = factor(species))
  stopifnot(all(species_targets$species %in% unique(rij_data$id)))
  rij_data <- rij_data %>% dplyr::select(-species,-area)
  
  feature_data <- dplyr::left_join(feature_data %>% mutate(id = factor(as.character(id))),
                                   species_targets %>% dplyr::rename(id = 'species'), by ='id' ) %>% 
    mutate(id = as.numeric(as.character(id)))

  print(paste0('No targets could be calculated for ', length( which(feature_data$id %notin% unique(rij_data$id)) ), ' species!'))
  
  # Correct name column and set to fname
  if(!assertthat::has_attr(feature_data,"name")) { feature_data$name <- NA}
  feature_data$name <- ifelse(is.na(feature_data$name),paste0(feature_data$id,"__",feature_data$binomial),feature_data$name)
  # Set columns for analysis
  rij_data$species <- rij_data$id; feature_data$species <- feature_data$id
  
  # Remove species not represented in PU table
  feature_data <- feature_data %>% drop_na(absolute_target)
  # Subset to relevant columns for the analysis and ignore rest
  feature_data <- feature_data %>% 
    dplyr::select(id,name,relative_target,species) %>% 
    distinct()
  # Reformat the ids and species to numbers
  rij_data$id <- as.numeric(as.character(rij_data$id))
  rij_data$species <- as.numeric(as.character(rij_data$species))
  rij_data$amount <- as.integer(rij_data$amount)
  gc()
}
myLog('Targets prepared!')
# Security checks
stopifnot( assert_that(all(unique(rij_data$species) %in% unique(feature_data$id) ))) # All species present in feature data?
#### Constrains and CO ####

if(ensure_all_land){
  # Include pu_data as feature
  rij_data <- bind_rows(rij_data,
                        pu_data %>% dplyr::select(id,cost) %>% dplyr::rename(pu = id, amount = cost) %>% 
                          dplyr::mutate(id = 77777777, species = 77777777))
  feature_data <- bind_rows(feature_data, data.frame(id = 77777777, species = 77777777,name = 'OVERALL.LAND',relative_target = 1,
                                                     stringsAsFactors = F))
}
# ------- #
# Load in Constraints data
constraints_data <- bind_rows(lapply(constraints_path, read_fst)) %>% 
  unique() %>% # Remove duplicated PU's
  dplyr::rename(protected = 'amount',id = 'pu') %>% 
  # Convert to 1000 range
  dplyr::mutate(protected = round( protected * 1000 ) ) %>% 
  # Work aound to ensure that constraints can't exceed cost
  left_join(., pu_data %>% dplyr::select(id,cost), by = "id" ) %>% 
  dplyr::mutate( protected = pmin(protected,cost) ) %>% 
  dplyr::select( - cost )

# Join back with PU data
pu_data <- dplyr::left_join(pu_data %>% mutate(id = factor(id)),
                            constraints_data %>% mutate(id = factor(id)),
                            by = 'id') %>% 
  dplyr::mutate(protected = ifelse(is.na(protected),0,protected) ) %>%  # If unprotected, assign 0
  dplyr::left_join(., regions_data %>% mutate(id = factor(id)), by = 'id') %>% # Also load in the biome/ecoregion ids
  dplyr::mutate(id = as.numeric(id)) # Convert back to numeric

# Append carbon data to rij-table and feature data    
if(include_carbon){
  # Rescale carbon to range of 1 to 1000
  carbon_data <- carbon_data %>% dplyr::filter(amount > 0)
  carbon_data$amount <- scales::rescale(carbon_data$amount,to = c(1,1000))
  # Join in Carbon data
  rij_data <- bind_rows(rij_data, carbon_data %>% dplyr::select(-id))
  feature_data <- bind_rows(feature_data,
                            tibble(name = "AGB_carbon", id = unique(carbon_data$id),
                                   species = ifelse(is.character(feature_data$species),"9999999",9999999),
                                   data = "AGB carbon", 
                                   size = n_distinct(carbon_data$pu),
                                   absolute_target = n_distinct(carbon_data$pu)*(target_resolution_number^2),relative_target = 1)
  )
  # Calculate carbon weights
  # Guarantee that carbon will be represented in the optimal solution 
  # by calculating the minimum weight necessary as the sum of the weight values 
  # for the other features
  if(phylo_weights) {stop('Phylo weights and carbon currently dont work!') }
  n = feature_data %>% dplyr::filter(id %notin% c(9999999,9999998)) %>% tally()
  cw <- (n$n+1) / carbon_multiplier
  carbonwater_weights <- c(rep(1,n$n), cw)
  rm(cw,n)
}
if(include_water){
  rij_data <- bind_rows(rij_data, water_data %>% dplyr::select(-id))
  feature_data <- bind_rows(feature_data,
                            tibble(name = "Water", id = unique(water_data$id),
                                   species = ifelse(is.character(feature_data$species),"9999998",9999998),                                   data = "Water", 
                                   size = n_distinct(water_data$pu),
                                   absolute_target = n_distinct(water_data$pu)*(target_resolution_number^2),relative_target = 1)
  )
  # Calculate water weights 
  # Guarantee that water will be represented in the optimal solution 
  # by calculating the minimum weight necessary as the sum of the weight values 
  # for the other features and adding a small number to it (e.g. 1).
  n = feature_data %>% dplyr::filter(id %notin% c(9999999,9999998)) %>% tally()
  ww <- (n$n+1) / water_multiplier
  if(include_carbon){ # Carbon also included?
    carbonwater_weights <- c(carbonwater_weights, ww)
  } else {
    carbonwater_weights <- c(rep(1,n$n), ww)
  }
  rm(ww,n)
}

# Include phylogenetic weights
if(phylo_weights){
  # Use the calculate phylo weights
  phylo_data <- read_fst(phylo_path) %>% dplyr::filter(iucn_id_no %in% feature_data$id)
  stopifnot( assert_that( all( has_name(phylo_data,c('ed','edge','rank_weight')) ),
                          any(duplicated(phylo_data$iucn_id_no)) == FALSE)
  )
  if(phylo_comparison){
    # If comparison should be made, take only those records which have estimates for both ED and TW
    phylo_data <- phylo_data %>% dplyr::filter(both_check == TRUE)
  }
  if(phylo_type == 'noweight'){
    phylo_data <- phylo_data[,c('iucn_id_no','rank_weight')] %>% drop_na() %>% distinct() %>% 
      # Reset all weight to 1 for a normal comparison
      dplyr::mutate(noweight = 1)
  } else {phylo_data <- phylo_data[,c('iucn_id_no',phylo_type)] %>% drop_na() %>% distinct()}

  # Limit feature ids only to those ids with a score
  feature_data <- feature_data %>% dplyr::filter(id %in% phylo_data$iucn_id_no)
  rij_data <- rij_data %>% dplyr::filter(id %in% phylo_data$iucn_id_no)
  assert_that(
    nrow(phylo_data) == nrow(feature_data),
    all( phylo_data$iucn_id_no %in% unique(rij_data$id) )
  )
  # Left_join with feature_data and set weighting score
  phylo_data <- phylo_data[,c('iucn_id_no',phylo_type)] %>% right_join(., feature_data,by = c('iucn_id_no' = 'id') )
  phylo_score <- phylo_data[,phylo_type]
  stopifnot( assert_that(length(phylo_score) == nrow(phylo_data)) )
}

# Minimimal Budget (e.g. already protected pu's)
min_budget <- (sum(pu_data$protected[which(pu_data$protected>0)]) / sum(pu_data$cost))
if(pa_lockedin) { myLog("Current amount of land protected: ", round(min_budget, 3) * 100, "%") }

# Some final checks that everything is in order
stopifnot(
  assert_that(
    nrow(pu_data) == pu_id_raster@data@max,
    all(rij_data$pu %in% pu_data$id),
    is.character(feature_data$name),
    anyNA(feature_data$name) == FALSE,
    anyNA(feature_data$id) == FALSE,
    assertthat::has_name(feature_data,"species"),
    anyDuplicated( feature_data$id )==FALSE,
    anyDuplicated( feature_data$name )==FALSE,
    assertthat::has_name(rij_data,"species"),
    is.numeric(rij_data$pu),
    !anyNA(rij_data$species),
    (round( sum(pu_data$protected),2 ) == round( min_budget * sum(pu_data$cost),2) )
  )
)

# Run variety if defined
if(run_variety== "carbonweights"){
  # Source the target script
  source('src/009_CarbonWeightProblems.R')
  stop('DONE')
} else if(run_variety== "carbonwaterweights"){
  if(include_carbon==FALSE || include_water == FALSE) stop('Carbon or water is not defined')
  # Source the target script
  source('src/009_CarbonWaterWeightProblems.R')
  stop('DONE')
} else if(run_variety == 'carbononly') {
  if(include_carbon==FALSE) stop('Carbon is not defined')
  # Source the target script
  source('src/009_CarbonOnlyProblems.R')
  stop('DONE')
}

myLog('Loaded overall data for ', n_distinct(rij_data$pu), ' PUs',' and ', n_distinct(feature_data$id), ' features')
# ------------------------------------------------------------------ #
set_number_of_threads( cores ) # Set parallel processing # probably obsolete
# Correct protected areas if set
if(pa_lockedin) {
  budgets <- sort( c(min_budget, budgets[which(budgets > min_budget)] ) ) 
  pu_data$cost <- pu_data$cost / 1000 ; pu_data$protected <- pu_data$protected / 1000
  rij_data$amount <- rij_data$amount / 1000 # Also set down the AOH by that amount
  if(ensure_all_land){
    rij_data$amount[rij_data$id==77777777] <- rij_data$amount[rij_data$id==77777777] / 1000 # Also set down the AOH by that amonut
  }
} else {budgets <- sort( budgets ) }

if(what=='minshort_speciestargets'){
#### Specific specific targets ####
  for(k in seq(1, length(budgets))){
    b = budgets[k]
    myLog('-----> Now processing budget option = ', round(b,3))
    
    # File name for adding carbon
    pa_carbon <- ifelse(include_carbon,"_carbon_","")
    # File name for adding water
    pa_water <- ifelse(include_water,"_water_","")
    # Filename for adding protected areas
    pa_fname <- ifelse(pa_lockedin,"withPA","")
    # Filename for splitid
    pa_split <- ifelse(split_species,paste0("_",split_id,"_"),"")
    # Is there a carbon weight multiplier greater than 1
    pa_carbmult <- ifelse(carbon_multiplier > 1, paste0("_carbweight",carbon_multiplier,"_"),"")
    # Is there a water weight multiplier greater than 1
    pa_watermult <- ifelse(water_multiplier > 1, paste0("_waterweight",water_multiplier,"_"),"")
    # Plants included?
    pa_plants <- ifelse(exclude_plants==TRUE,"_waPlants","")
    # Phylogenetic weight
    pa_phylo <- ifelse( phylo_weights == TRUE, paste0("_phylo-",phylo_type,ifelse( phylo_comparison == TRUE,'-comparison',''),"_"), "")
    
    # Addition if this is a representative subsample
    if(repr_id){
      pa_repr_add <- paste0('REPrun_',permut_run,"/")
      dir.create(paste0(results_path,'/',pa_repr_add),showWarnings = FALSE)
    } else { pa_repr_add <- "" }
    
    ## Construct output file names
    # .Platform$file.sep # For a system specific file separator
    out_name <- paste0(results_path,'/',pa_repr_add,'minshort_speciestargets',pa_split,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,'_',round(budgets[k],2)*100,'perc','.fst')
    out_name_tif <- paste0(results_path,'/',pa_repr_add,'minshort_speciestargets',pa_split,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,'_',round(budgets[k],2)*100,'perc','.tif')
    out_name_save <- paste0(results_path,'/',pa_repr_add,'minshort_speciestargets_securitysave',pa_split,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,'_',round(budgets[k],2)*100,'perc','.fst')
    out_previous_save <- paste0(results_path,'/',pa_repr_add,'minshort_speciestargets_securitysave',pa_split,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,'_',round(budgets[k-1],2)*100,'perc','.fst')
    
    # -------------------- #
    if((!overwrite_existing_data) && file.exists(out_name)){
      print('Prioritization results do already exist -> Skip')
      next()
    } 
    # -------------------- #
    # Lock in previous one if higher than 1
    if(k>1){
      assert_that(file.exists(out_previous_save));rm(prev_solution)
      prev_solution <- read_fst(out_previous_save)
      prev_solution <- subset(prev_solution,select = c("id","solution_1")) %>% dplyr::rename(pu = "id",lower = "solution_1") %>% 
        dplyr::filter(lower > 0) %>% # Filter out unselected PUs
        dplyr::mutate(upper = 1) %>% # Add upper bound of whatever the maximum is 
        dplyr::mutate(lower = pmin(lower, 1 )) # Assert that lower does not overshoot 1
      # Lock in previous solution if it was binary or as continious
      stopifnot( assert_that(nrow(prev_solution)>0 ) )
    }
    
    # Budget formulation
    b_cells <- b * sum(pu_data$cost)
    # create problem
    start_time <- Sys.time()
    p_maxfeat <- problem(pu_data, feature_data, rij_data, cost_column = "cost") %>%
      add_min_shortfall_objective(budget = b_cells) %>% 
      add_relative_targets('relative_target') %>%
      add_gurobi_solver(gap = optimality_gap,time_limit = time_limit,threads = cores,numeric_focus = FALSE,
                        verbose = TRUE) 
    # Decision type
    if(decision_type=="binary") {
      p_maxfeat <- p_maxfeat %>% add_binary_decisions()
    } else{
      p_maxfeat <- p_maxfeat %>% add_proportion_decisions()
    }
    # Locked in constrains and feature weights
    if(k>1) { p_maxfeat <- p_maxfeat %>% add_manual_bounded_constraints( prev_solution ) } # Manually lock in previous solution
    if(pa_lockedin & k == 1){ 
      p_maxfeat <- p_maxfeat %>% add_manual_locked_constraints(
        pu_data %>% dplyr::filter(protected > 0) %>% 
          dplyr::select(id,protected) %>%
          dplyr::rename(pu = "id", status = "protected")
      )
    }
    # Include phylo
    if(phylo_weights){ p_maxfeat <- p_maxfeat %>% add_feature_weights(phylo_score) }
    if(include_carbon | include_water ){ p_maxfeat <- p_maxfeat %>% add_feature_weights(carbonwater_weights)  } # Add carbon &/ water
    # Add portfolio to the solution
    if(port){ p_maxfeat <- p_maxfeat %>% add_pool_portfolio(method = 0) } # Using method 2 allows to obtain solutions close to optimality
    
    # Solving   
    r_maxfeat <- prioritizr::solve(p_maxfeat, force = F,run_checks = FALSE)
    end_time <- Sys.time(); print( round((end_time - start_time),2) )
    # --- #
    myLog('Assess feature representation..')
    # Seucrity save 
    write_fst(r_maxfeat,out_name_save) 
    
    # Calculate representation for top solution
    if(calculate_representation){
      out <- feature_representation2(p_maxfeat,r_maxfeat[, "solution_1", drop = FALSE], cores_repr)
      # Join in amount necessary for reaching target as well as feature abundance in planning_units()
      out$absolute_target <- p_maxfeat$targets$output()$value
      out$feature_abundance_pu <- p_maxfeat$feature_abundances_in_planning_units()
      # Append target
      write.fst(out, out_name )
    }
    # Create raster output
    out_ras <- createOutput(pu_id_raster,r_maxfeat,"solution_1")
    names(out_ras) <- paste0("percent",round(b,2) * 100) # Rename
    # Save output
    writeGeoTiff(out_ras, out_name_tif,dt = ifelse(decision_type=="binary","INT2S","FLT4S") )
    rm(r_maxfeat,p_maxfeat,out);gc()
  } # End of budgets calculation
  rm(carbonwater_weights)
} # Minshort statement end

}

stop("DONE!")
