# ------------------------------------------------------- #
# Preparation script for all IUCN data
library(raster)
library(fasterize)
library(rgdal)
library(velox)
library(sf)
library(stringr)
library(assertthat)
library(doParallel)
source("src/000_ConvenienceFunctions.R")
# Adjust temporary folder
tmp <- '/media/martin/data/temp'
dir.create(tmp,showWarnings = FALSE)
rasterOptions(tmpdir = tmp)
Sys.setenv('TMP' = '/media/martin/data/temp')

use_velox <- TRUE # use the Velox package
# Number of cores and cluster_type for parallel processing 
cores <- 4

# The Path to the biodiversity data
biodiversity_list = "/media/martin/data/raw/IUCN_RL_2019_2_SpeciesList.csv"
final_path <- "/media/martin/data/features/"

# Global grid path
gridpath = "data/"
projection = "mollweide"
resolution = "10km"
target_resolution_number <- as.numeric(gsub("\\D","",resolution)) # Target resolution as number 

# Classes #
taxa <- c('AMPHIBIA','AVES','MAMMALIA','REPTILIA','PLANTAE')

species <- readr::read_csv(biodiversity_list)

stopifnot(
  assert_that(
    is.numeric(cores),
    file.exists(biodiversity_list),
    dir.exists(gridpath),
    dir.exists(final_path),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif")),
    # The biodiversity files
    file.exists( paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",taxa[1],".gpkg") ),
    file.exists( paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",taxa[2],".gpkg") ),
    file.exists( paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",taxa[3],".gpkg") ),
    file.exists( paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",taxa[4],".gpkg") ),
    file.exists( paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",taxa[5],".gpkg") ),
    # Check that target groups are present
    all( taxa[-5] %in% species$class ),
    all( taxa[5] %in% species$kingdom )
  )
)

# ------------------------------------------------------- #
# Make a global land area map out of the grid
globalgrid <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif"))
globalgrid1km <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif"))
globalgrid1km[globalgrid1km>=1] <- 1 # Make a global mask
writeGeoTiff(globalgrid1km,paste0(tmp,"/globalgrid1km.tif"))
globalgrid1km <- raster(paste0(tmp,"/globalgrid1km.tif"))

# And the fraction cover
globalgrid_fraction <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif")) / 1000
writeGeoTiff(globalgrid_fraction,paste0(tmp,"/globalgrid_frac.tif"),'FLT4S')
globalgrid_fraction <- raster(paste0(tmp,"/globalgrid_frac.tif"))

stopifnot(
  assert_that(
    res(globalgrid1km)[1] == (target_resolution_number*1000),
    res(globalgrid_fraction)[1] == (target_resolution_number*1000)
  )
)

if(use_velox){ vx <- velox(globalgrid1km) } else{ vx <- NULL} # Convert to velox

# Fire up cores
if(file.exists('log.txt')) file.remove('log.txt')
cl <- makePSOCKcluster(cores, outfile = "log.txt")
registerDoParallel(cl)

# -------------------------------------- #
# Now loop over all groups
for(target_taxa in taxa) {
  myLog('---------------------------')
  myLog('|-> Processing ',target_taxa, ' at ',resolution)
  
  # Input and Output path
  biodiversity_path = paste0("/media/martin/data/raw/IUCN_RL_2019_2_Species_Cleaned_",target_taxa,".gpkg")
  output_path = paste0(final_path,target_taxa,"_",resolution)
  # Some error checks
  if(!dir.exists(output_path)) { dir.create(output_path)}
  
  # # Load in the biodiversity data
  biodiversity <- st_read(biodiversity_path,
                          stringsAsFactors = FALSE,
                          promote_to_multi = TRUE)
  
  # --- Subsetting --- #  
  if(target_taxa == "PLANTAE"){
    biodiversity <- biodiversity %>% dplyr::filter(kingdom == target_taxa)
  } else {
    biodiversity <- biodiversity %>% dplyr::filter(class == target_taxa)
  }
  # Filter metadata
  biodiversity <- biodiversity %>% 
    # Remove extinct species
    dplyr::filter(category %notin% c("EX","EW")) %>% 
    # Only species with 1 (Extant) and 4 (possibly Extinct)
    dplyr::filter(presence %in% c(1, 4)) %>% 
    # Only species with native or reintroduced origin
    dplyr::filter(origin %in% c(1,2)) %>% 
    # Only species being resident (1), Breeding (2), Non-breedin season (3), passage (4) or (5) Seasonal Occurrence Uncertain
    dplyr::filter(seasonal %in% c(1,2,3,4,5)) %>% 
    # Make fname
    #  dplyr::left_join(., species %>% dplyr::select(id_no,family), by = 'id_no') %>% 
    dplyr::mutate(fname = paste0(str_to_lower(family),'_',str_replace_all(str_to_lower(binomial)," ","_") ) )
  # If Plantae use binomial only to allign with other plant data
  if(target_taxa == "PLANTAE"){ biodiversity$fname <- str_replace_all(str_to_lower(biodiversity$binomial)," ","_") }
  
  # Reproject all
  if(projection!= "wgs84"){
    biodiversity <- st_transform(biodiversity,crs = proj4string(globalgrid))
  }
  # Now loop over all files
  o <- foreach::foreach( species = iter(unique(biodiversity$fname) ),
                         .inorder = T,
                         .errorhandling = 'pass',
                         .packages = c("fasterize","raster","stringr","velox","assertthat"),
                         .export = c("biodiversity","globalgrid1km","globalgrid","writeGeoTiff","resolution","alignRasters","output_path",
                                     "vx","use_velox","globalgrid_fraction")
  ) %dopar% {
    sub <- subset(biodiversity,fname == species)
    fname = species
    if(file.exists(paste0(output_path,"/",fname,".tif"))) { return(NULL) }
    
    myLog("Rasterizing at ",str_to_title(species))
    if(use_velox){
      # Use the velox raster package
      sub$dummy <- 1
      vx2 <- vx$copy()
      ## Rasterize, set background to -1
      try( vx2$rasterize(spdf=sub, field="dummy", background=0,small = TRUE), silent = TRUE)
      species_range <- vx2$as.RasterLayer()
      rm(vx2)
    } 
    if(!exists('species_range') | use_velox == FALSE) {
      # Use the fasterize raster package
      species_range <- fasterize::fasterize(sf = sub,
                                 raster = globalgrid1km,
                                 fun = "any",
                                 field = NULL, # Sets 1 to all cells covered
                                 background = 0 # Set to zero for aggregation!
      )
      # Now aggregate to target resolution
      species_range <- raster::aggregate(species_range, fact = target_resolution_number, na.rm = TRUE, fun = mean)
    }
    # And align
    if(!compareRaster(species_range,globalgrid_fraction)){
      species_range <- alignRasters(species_range,globalgrid,method = 'ngb',func = mean,cl = FALSE)
    }
    # Remove fractional cover
    species_range <- (species_range * globalgrid_fraction)
    species_range[species_range == 0] <- NA
    # Clamp and convert to Integer
    species_range <- raster::clamp(species_range, lower = 0, upper = 1)
    species_range <- round( species_range * 1000 )
    
    # Save outputs 
    writeGeoTiff(species_range,
                 fname = paste0(output_path,"/",fname),
                 dt = "INT2S"
    ) 
    # Clean up files
    rm(species_range)
    removeTmpFiles(.5)
    #ll <- list.files(tempdir(),pattern = "*.gri",full.names = TRUE)
    #unlink(ll);rm(ll)
  }
  # Save output list to check later
  saveRDS(o, paste0(final_path,target_taxa,"_",resolution,".rds") )
  
  gc()
}

stopCluster(cl)
stop("DONE!")

assert_that(all(is.null( unlist(o)) ))
