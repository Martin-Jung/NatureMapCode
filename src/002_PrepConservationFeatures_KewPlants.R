# ------------------------------------------------------- #
# The preparation script for the Kew digitized plant data
# Creates rasterized species range layers per species
# Using a template grid as target
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
biodiversity_path = c('/media/martin/data/raw/Plantdata_Kew/')
biodiversity_list = "/media/martin/data/raw/IUCN_RL_2019_2_SpeciesList.csv"

# Global grid path
gridpath = "data/"
projection = "mollweide"
resolution = "10km"
target_resolution_number = 10

# Dataset
ds <- c('IUCN','BGCI')

# Some error checks
stopifnot(
  assert_that(
    all(dir.exists(biodiversity_path)),
    file.exists(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif"))
  )
)

# ------------------------------------------------------- #
# Make a global land area map out of the grid
globalgrid <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif"))
globalgrid[globalgrid>=1] <- 1 # Make a global mask
writeGeoTiff(globalgrid,paste0(tmp,"/globalgrid.tif"))
globalgrid <- raster(paste0(tmp,"/globalgrid.tif"))

# And the fraction cover
globalgrid_fraction <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif")) / 1000
writeGeoTiff(globalgrid_fraction,paste0(tmp,"/globalgrid_frac.tif"),'FLT4S')
globalgrid_fraction <- raster(paste0(tmp,"/globalgrid_frac.tif"))

if(use_velox){ vx <- velox(globalgrid) }
# Load in the biodiversity data
species <- readr::read_csv(biodiversity_list) %>% 
  # Immediate filter
  dplyr::filter(kingdom == 'PLANTAE') %>% 
  dplyr::filter(category %notin% c("EX","EW")) # Remove extinct species

# Fire up cores 
if(file.exists('log.txt')) file.remove('log.txt')
cl <- makePSOCKcluster(cores, outfile = "log.txt")
registerDoParallel(cl)

# Loop over sources
for(src in ds ){
  myLog('---------------------------')
  myLog('|-> Processing ',src, ' at ',resolution)
  
  # List of files
  if(src == "BGCI"){
    biodiversity <- c(
      list.files(paste0(biodiversity_path,src,"/polygons"),pattern = "*.shp",recursive = T,full.names = T),
      list.files(paste0(biodiversity_path,src,"/NEW/Outputs/polygons"),pattern = "*.shp",recursive = T,full.names = T)
    )
  } else {
    biodiversity <- c(
      list.files(paste0(biodiversity_path,src,"/polygons"),pattern = "*.shp",recursive = T,full.names = T)
    )
  }
  # Remove species not in the IUCN species list
  biodiversity <- biodiversity[which( tools::file_path_sans_ext(basename(biodiversity)) %in% species$binomial )]
  
  # Output path
  output_path = paste0("/media/martin/data/features/KEW_PLANTAE_",src,"_",resolution)
  if(!dir.exists(output_path)){ dir.create(output_path) }
  
  # --- Processing start ---- #
  o <- foreach::foreach( species = iter(biodiversity),
                         .inorder = F,
                         .packages = c("fasterize","raster","stringr","velox","assertthat","sf"),
                         .export = c("biodiversity","projection","globalgrid","globalgrid_fraction","writeGeoTiff","vx","use_velox")
  ) %dopar% {
    #for(species in unique(biodiversity$fname)) {
    sub <- sf::st_read(species,quiet = T) 
    fname = tools::file_path_sans_ext(basename(species))
    if(file.exists(paste0(output_path,"/",fname,".tif"))) { return(NULL) }
    
    if(projection!= "wgs84"){
      sub <- st_transform(sub,crs = proj4string(globalgrid))
    }
    sub$dummy <- 1 # Make a dummy field
    
    if(use_velox){
      # Use the velox raster package
      vx2 <- vx$copy()
      ## Rasterize, set background to -1
      vx2$rasterize(spdf=sub, field="dummy", background=NA,small = TRUE)
      species_range <- vx2$as.RasterLayer()
      rm(vx2)
      
    } else {
      # Use the fasterize raster package
      species_range <- fasterize(sf = sub,
                                 raster = globalgrid,
                                 fun = "any",
                                 field = NULL, # Sets 1 to all cells covered
                                 background = NA
      ) # Fasterize
    }
    assert_that(
      cellStats(species_range,'max') > 0
    )
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
    rm(species_range)
    return(NULL)
  }
  # Save output list to check later
  saveRDS(o, paste0('/media/martin/data/features/',"KEW_",src,"_",resolution,".rds") )
  
}
stopCluster(cl);stopImplicitCluster()
stop("DONE!")
