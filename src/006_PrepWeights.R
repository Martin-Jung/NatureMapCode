# Purpose of this script is to prepare weights for the prioritization 
# analysis
set.seed(3263)

# IUCN threat status and categories list
biodiversity_list = "/media/martin/data/raw/IUCN_RL_2019_2_SpeciesList.csv"
output_path = "data"

# Paths to phylogenetic data
phylopath_birds = "/media/martin/data/weights_phylogeny/birds/"
phylopath_amphibians = "/media/martin/data/weights_phylogeny/amphibians/"
phylopath_squamates = "/media/martin/data/weights_phylogeny/squamates/"
phylopath_mammals = "/media/martin/data/weights_phylogeny/mammals/"

# EDGE scores
edge_scores_path = "/media/martin/data/weights_phylogeny/EDGE_Lists_2019.xlsx"
# Plant threat status
plant_threat <- "/media/martin/data/weights_threat/ThreatenedPlants140120.xlsx"
  
# IUCN API key (doesn't get pushed to github ;)  )
iucn_api_key = suppressWarnings( read.table("data/key.txt",stringsAsFactors = F)[,1] )

# --------- Packages ------------- #
library(tidyverse)
library(assertthat)
library(stringr)
library(fst)
library(ape)
library(picante)
library(phytools)
library(phylobase)
library(taxize)
library(geiger)
library(readxl)
library(rredlist)
library(rgbif)
source("src/000_ConvenienceFunctions.R")

# Get Species list. Names in column "binomial"
species <- read_fst("data/speciesID_table_10km.fst") %>% dplyr::select(-fname)
species$kingdom[species$data=='GARD Shai'] <- 'ANIMALIA'
species <- species #%>% dplyr::filter(data %in% c('GARD Shai','IUCN Amphibians','IUCN Bird',
                   #                             'IUCN Mammal','IUCN Reptiles'))
# --------------------------------------------------------- #
#### Format EDGE scores for all species ####
# Load ED scores
ed_amph <- readxl::read_xlsx(edge_scores_path,sheet = "ED scores - amphibians") %>% 
  dplyr::select(1:6) %>% dplyr::rename(iucn_id_no = `RL ID`, category = `RL Status`)
ed_bird <- readxl::read_xlsx(edge_scores_path,sheet = "ED scores - birds" ) %>% 
  dplyr::select(1:4,6:7) %>% dplyr::rename(iucn_id_no = `RL ID`, category = `RL Status`) %>% mutate(iucn_id_no = as.character(iucn_id_no))
ed_mam <- readxl::read_xlsx(edge_scores_path,sheet = "ED scores - mammals" ) %>% 
  dplyr::select(1:6) %>% dplyr::rename(iucn_id_no = `RL ID`, category = `RL Status`) %>% mutate(iucn_id_no = as.character(iucn_id_no))
ed_rep <- readxl::read_xlsx(edge_scores_path,sheet = "ED scores - reptiles" ) %>% 
  dplyr::select(1:6) %>% dplyr::rename(iucn_id_no = `RL ID`, category = `RL Status`)
# Correct names
names(ed_amph) <- str_to_lower(names(ed_amph))
names(ed_bird) <- str_to_lower(names(ed_bird))
names(ed_mam) <- str_to_lower(names(ed_mam))
names(ed_rep) <- str_to_lower(names(ed_rep))

# Combined and format all names
ed2019 <- bind_rows(
  ed_amph %>% mutate(class = 'AMPHIBIA'),
  ed_bird %>% mutate(class = 'AVES'),
  ed_mam %>% mutate(class = 'MAMMALIA'),
  ed_rep %>% mutate(class = 'REPTILIA')
) %>% 
  dplyr::mutate(iucn_id_no = as.numeric(iucn_id_no) ) %>% 
  dplyr::select(-category) %>% 
  drop_na(species) %>% 
  distinct() # Some species are duplicate in here :-(
rm(ed_amph,ed_bird,ed_mam,ed_rep) # Clean up
# Security check
stopifnot(assert_that(
  any(duplicated(ed2019$species))==FALSE
))
# ---- # 
# Match id's with iucn_id list
ed2019$id <- 1:nrow(ed2019)
# Query the species not already in species list
ed_match <- ed2019 %>% dplyr::filter(!is.na(iucn_id_no)) %>% 
  inner_join(., species %>% dplyr::select(iucn_id_no), by = "iucn_id_no")

# All those still missing (e.g. NonMatched binomial, NA)
miss <- ed2019[which(ed2019$id %notin% ed_match$id),] 
assert_that( nrow(miss) + nrow(ed_match) == nrow(ed2019))

myLog('Remaining number of unmatched species: ', nrow(miss))
# ----------------------------------------- #
# 1) 
# Check which species are already in species list with binomial name
nonmiss <- miss[which( str_to_lower(miss$species) %in% str_to_lower(species$binomial) ),]
miss <- miss[which( str_to_lower(miss$species) %notin% str_to_lower(species$binomial) ),] # The remaining species 

# Fill in / overwrite ids with those found in the species_id list
nonmiss$iucn_id_no <- species$iucn_id_no[match(str_to_lower(nonmiss$species),str_to_lower(species$binomial))]
# Join with match and clean up
ed_match <- bind_rows(ed_match,nonmiss)
assert_that( nrow(miss) + nrow(ed_match) == nrow(ed2019),
             any(duplicated(ed_match$iucn_id_no))==FALSE)

myLog('Remaining number of unmatched species: ', nrow(miss))

# Also get the species not in ed_match and miss
species_nonmatch <- species %>% dplyr::filter(str_to_lower(binomial) %notin% unique( c( str_to_lower(ed_match$species), str_to_lower(nonmiss$species) ) ) )
species_nonmatch <- species_nonmatch %>% dplyr::filter(data %in% c('GARD Shai','IUCN Amphibians','IUCN Bird',
                                                                        'IUCN Mammal','IUCN Reptiles')) %>% 
  dplyr::select(iucn_id_no,data,order,family,genus,binomial,class)
write_csv(species_nonmatch,'ED_notpresent.csv')
rm(nonmiss);gc()
# ----------------------------------------- #
# 2)
# Match remaining missing names against gnr
pb <- progress::progress_bar$new(total = nrow(miss))
synonym1 <-  data.frame( species = character() ,synonyms = list(),stringsAsFactors = FALSE)
# Loop through missing matches
for(r in 1:nrow(miss)){
  #n <- taxize::gnr_resolve(miss$species[r],resolve_once = F, with_context = FALSE)
  # 'Ambystoma subsalsum' --> Ambystoma velasci 
  # First query IUCN 
  n <- synonyms_df(
    taxize::synonyms(miss$species[r],db='iucn',key=iucn_api_key)
  )
  if(nrow(n)>0 & has_name(n,'syn_name')){
    synonym1 <- bind_rows(
      synonym1,
      data.frame(species = miss$species[r], syn_name = n$syn_name,stringsAsFactors = FALSE)
    )
  } else {
    n <- synonyms_df(
      taxize::synonyms(miss$species[r],db='itis')
    )
    if(nrow(n)>0 & has_name(n,'syn_name')){
      synonym1 <- bind_rows(
        synonym1,
        data.frame(species = miss$species[r], syn_name = n$syn_name,stringsAsFactors = FALSE)
      )
    }
  }
  pb$tick()
}
myLog('Found synonyms for ', n_distinct(synonym1$species),' out of ', n_distinct(miss$species))
#saveRDS(synonym1,'synonym.rds')
synonym1 <- readRDS('synonym.rds')
# Get Species matches
miss_species <- species %>% dplyr::filter(iucn_id_no %notin% ed2019$iucn_id_no)

# Split synonyms into binomial pairs to ensure species
synonym1 <- synonym1 %>% tidyr::separate(syn_name,into = c('genus','syn_species','ssp'),sep = ' ') %>% 
  dplyr::mutate(syn_name_clean = paste0(genus," ", syn_species)) # Combine into a species synonym again
# Join in ED
synonym1 <- left_join(synonym1,miss %>% dplyr::select(species,ed), by = c('syn_name_clean' = 'species') )
synonym1 <- synonym1 %>% dplyr::select(-genus,-syn_species,-ssp) %>% distinct() %>% drop_na()

# Make some checks
any( synonym1$species %in% species$binomial )
any( synonym1$syn_name_clean %in% species$binomial )
# Both not present...

myLog('Remaining number of unmatched species: ', nrow(miss))
rm(syn_match)
# ----------------------------------------- #
# Correct those two missing ones
ed2019[which(ed2019$species == 'Atopophrynus syntomopus'),c('order','family')] <- c('ANURA', 'CRAUGASTORIDAE')
ed2019[which(ed2019$species == 'Geobatrachus walkeri'),c('order','family')] <- c('ANURA', 'CRAUGASTORIDAE')
assert_that(anyNA(ed2019$family)==FALSE,anyNA(ed2019$genus)==FALSE)

# 3)
# Calculate genus averages for each genus in the full dataset
ed_avggenus <- ed2019 %>% tidyr::separate(species,into = c('genus','species','ssp'),sep = ' ') %>% 
  dplyr::mutate(genus = str_to_title(genus)) %>% 
  dplyr::group_by(family,genus) %>% 
  dplyr::summarise(ed = mean(ed,na.rm = T))
# And per family
ed_avgfamily <- ed2019 %>% 
  dplyr::mutate(family = str_to_title(family)) %>% 
  dplyr::group_by(family) %>% 
  dplyr::summarise(ed_fam = mean(ed,na.rm = T))

# Now match with species average
ed_avgmatch <- species_nonmatch %>% 
  dplyr::filter(data %in% c('GARD Shai','IUCN Amphibians','IUCN Bird',
                                                         'IUCN Mammal','IUCN Reptiles')) %>% 
  dplyr::mutate(genus = str_to_title(genus), family = str_to_title(family)) %>% 
  dplyr::left_join(., ed_avggenus, by = c('family','genus')) %>% 
  # And family level
  dplyr::left_join(., ed_avgfamily, by = 'family') %>% 
  # If genus not existing, use
  dplyr::mutate(ed = ifelse(is.na(ed),ed_fam,ed)) %>% dplyr::select(-ed_fam) %>% 
  dplyr::mutate(genusfamavg = TRUE)
rm(ed_avggenus,ed_avgfamily)

# Remove all from average that are already in ed_match
ed_avgmatch <- ed_avgmatch %>% dplyr::filter(iucn_id_no %notin% ed_match$iucn_id_no)

# Merge the family_genus averaged ones into here 
ed_out <- bind_rows(
  ed_match %>% dplyr::mutate(genusfamavg = FALSE) %>% dplyr::rename(binomial = species) %>% dplyr::select(iucn_id_no, binomial, ed,genusfamavg),
  ed_avgmatch %>% dplyr::select(iucn_id_no, binomial, ed,genusfamavg)
) %>% drop_na()
assert_that( any( duplicated( ed_out$iucn_id_no ) )== FALSE )

# How many filled?
table(ed_out$genusfamavg)

# ----------------- #
# Following Montesino Pouzols et al. (2014)
tw <- data.frame(category = c("CR","EN","VU","NT","LC","LR/cd","LR/lc","LR/nt","DD",'Threatened'),
                 rank_weight = c(8,6,4,2,1,1,1,2,2,6))

# Join in categories and reformat them to numeric
ed_out <- ed_out %>% left_join(., species %>% 
                                 dplyr::select(iucn_id_no,category) %>% 
                                 left_join(.,tw, by = "category"), 
                               by = 'iucn_id_no')
assert_that( any( duplicated( ed_out$iucn_id_no ) )== FALSE )
# --- # 
# Join in all other species with threat data
# First with existing threat status
ed_rest1 <- species %>% dplyr::filter(iucn_id_no %notin% ed_out$iucn_id_no) %>% 
  dplyr::select(iucn_id_no,binomial,category) %>% drop_na() %>% distinct() %>% 
  left_join(.,tw, by = "category")
assert_that(any(is.na(ed_rest1$rank_weight))==FALSE,
            any(duplicated( ed_rest1$iucn_id_no )==FALSE))

# Now load in the BGCI plant threat data
tw_plant <- read_xlsx( plant_threat,sheet = 2 ) %>% 
  dplyr::rename(iucn_id_no = ID, binomial = TaxonName, category = ConsAssCategory) %>% 
  dplyr::filter(iucn_id_no %notin% ed_rest1$iucn_id_no)

# First those with id
ed_rest2 <- tw_plant %>% dplyr::filter(iucn_id_no %notin% c(ed_rest1$iucn_id_no)) %>% 
  dplyr::select(iucn_id_no, binomial, category) %>% 
  drop_na() %>% 
  left_join(.,tw, by = "category")
assert_that( any( duplicated( ed_rest2$iucn_id_no ) )== FALSE )

# Then match the remaining names to the binomial name where this is possible
ed_rest3 <- tw_plant %>% dplyr::filter(iucn_id_no %notin% c(ed_rest1$iucn_id_no,ed_rest2$iucn_id_no) ) %>% 
  dplyr::select(binomial, category) %>% 
  drop_na() %>% 
  # left join with species based on binomial
  left_join(., species %>% dplyr::select(iucn_id_no, binomial), by = 'binomial') %>% drop_na() %>% 
  left_join(.,tw, by = "category")
assert_that( any( duplicated( ed_rest3$iucn_id_no ) )== FALSE )

# Combine them all
ed_rest <- bind_rows(
  ed_rest1, ed_rest2, ed_rest3
)
rm(ed_rest1,ed_rest2,ed_rest3)
assert_that( any( duplicated( ed_rest$iucn_id_no ) )== FALSE )
ed_out <- bind_rows(
 ed_out, ed_rest %>% dplyr::filter(iucn_id_no %notin% ed_out$iucn_id_no)
)
assert_that( any( duplicated( ed_out$iucn_id_no ) )== FALSE )

# And finally all species without categories
ed_rest4 <- species %>% dplyr::filter(binomial %notin% ed_out$binomial) %>% 
  dplyr::select(iucn_id_no, binomial, category) %>% distinct()
ed_rest4$category <- 'DD';ed_rest4$rank_weight <- 2
assert_that( all(!is.na( ed_rest4$category )), any( duplicated( ed_rest4$iucn_id_no ) ) == FALSE )  

# Join in
ed_out <- bind_rows(ed_out, ed_rest4 %>% dplyr::filter(binomial %notin% ed_out$binomial, iucn_id_no %notin% ed_out$iucn_id_no) )

# FIX and the end
ed_out <- ed_out[which(ed_out$iucn_id_no %in% species$iucn_id_no),]

assert_that( anyNA(ed_out$iucn_id_no) == FALSE,
             anyDuplicated(ed_out$iucn_id_no)==0,
             all(ed_out$iucn_id_no %in% species$iucn_id_no),
             nrow(ed_out) <= nrow(species))

# Calculate EDGE score
ed_out <- ed_out %>% dplyr::mutate(edge = log((1+ed) * (2^rank_weight) ) )

# Create a check in the data.frame for species with both TW and ED
ed_out$both_check <- !is.na(ed_out$edge)

myLog('Number species with ED data: ',length(ed_out$ed[which(!is.na(ed_out$ed))]),
      '| Number species with threat scores: ', length(ed_out$rank_weight[which(!is.na(ed_out$rank_weight))]) )

write_fst(ed_out, "data/ED_threat_scores.fst")

rm(ed,ed2)
stop("DONE!")
# --------------------------------------------------------- #