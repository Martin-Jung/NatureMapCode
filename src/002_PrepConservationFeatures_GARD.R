# ------------------------------------------------------- #
# The preparation script for the GARD data 
# Creates rasterized species range layers per species
# Using a template grid as target
library(raster)
library(fasterize)
library(rgdal)
library(velox)
library(taxize)
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

cores <- 4 # Number of cores to use
use_velox <- TRUE # use the Velox package

# The Path to the biodiversity data
biodiversity_path = "/media/martin/data/raw/reptiles_gard/GARD1.5/gard_1_5_ranges_Maria.shp"
id_path = "data/newGardREP_IDs.csv" # List from Andy
# Global grid path
gridpath = "data/"
projection = "mollweide"
resolution = "10km"
target_resolution_number = 10
exportShape = FALSE # Should the shapefile be exported?

# Output path
output_path = paste0("/media/martin/data/features/GARDreptilia_",resolution)

# Some error checks
if(!dir.exists(output_path)) dir.create(output_path)
assert_that(
  file.exists(biodiversity_path),
  file.exists(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif")),
  dir.exists(output_path)
)

# ------------------------------------------------------- #
# Make a global land area map out of the grid
globalgrid <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif"))
globalgrid <- setMinMax(globalgrid)
globalgrid[globalgrid>=1] <- 1
writeGeoTiff(globalgrid,paste0(tmp,"/globalgrid.tif"))
globalgrid <- raster(paste0(tmp,"/globalgrid.tif"))

# And the fraction cover
globalgrid_fraction <- raster(paste0(gridpath,"globalgrid_",projection,"_",resolution,".tif")) / 1000
writeGeoTiff(globalgrid_fraction,paste0(tmp,"/globalgrid_frac.tif"),'FLT4S')
globalgrid_fraction <- raster(paste0(tmp,"/globalgrid_frac.tif"))

if(use_velox){ vx <- velox(globalgrid) }
# Load in the biodiversity data
biodiversity <- st_read(biodiversity_path)
biodiversity$fname <- str_replace_all(str_to_lower(biodiversity$BINOMIAL)," ","_")
# Export GARD taxon id's
ids <- read.csv(id_path,stringsAsFactors = FALSE)
biodiversity <- dplyr::left_join(biodiversity, ids %>% dplyr::select(gard_id,binomial), by = c('BINOMIAL' = 'binomial'))
assert_that(anyNA(biodiversity$gard_id)==FALSE)
#biodiversity$gard_id = as.numeric(str_remove(biodiversity$TAXONID,"R")) + 1000000000

if(exportShape){
  biodiversity <- biodiversity %>% dplyr::select(TAXONID,BINOMIAL,GROUP,gard_id) %>% 
    dplyr::rename(id_no = TAXONID,group = GROUP,binomial = BINOMIAL) %>%  dplyr::mutate(data = "gard_reptiles")
  biodiversity <- lwgeom::st_make_valid(biodiversity) # Make valid
  write_sf(biodiversity,"/media/martin/data/gard_reptiles.shp")
  
  # Get ID's
  # match all ed names against gnr
  n <- classification( as.character( biodiversity$binomial ), db = 'gbif',ask = F)
  matched_taxonomy <- data.frame()
  for(s in 1:length(n)){
    if(!is.null( nrow( n[[s]] ) )){ 
        sub <- n[[s]] %>% dplyr::select(name,rank) %>% dplyr::filter(rank %in% c('class','order','family')) %>% 
          tidyr::spread(rank,name) %>% 
          dplyr::mutate( binomial = names(n[s]) )
        
    } else {
      sub <- data.frame(class = 'Reptilia', family = NA, order = NA, binomial = names(n[s])) 
    }
      matched_taxonomy <- dplyr::bind_rows(
        matched_taxonomy,
        sub
      )
  }
  assert_that( anyNA(matched_taxonomy$binomial)==FALSE )
  # That one did not work for some reason
  #matched_taxonomy <- dplyr::bind_rows(
  #  matched_taxonomy,
  #  data.frame(class = 'Reptilia', family = NA, order = NA, binomial =  'Wallaceophis gujaratensis' )
  #)
  assert_that(nrow(matched_taxonomy) == nrow(biodiversity),
              any(!duplicated(matched_taxonomy$binomial)))
  
  x <- as.data.frame(biodiversity) %>% dplyr::select(-geometry)
  x <- dplyr::left_join(x, matched_taxonomy, by = 'binomial')
  saveRDS(x,'data/gard_ids.rds')
}

if(projection!= "wgs84"){
  biodiversity <- st_transform(biodiversity,crs = proj4string(globalgrid))
}

# Fire up cores
cl <- makePSOCKcluster(cores, outfile = "log.txt")
registerDoParallel(cl)

o <- foreach::foreach( species = iter(unique(biodiversity$fname) ),
                  .errorhandling = 'pass',
                  .packages = c("fasterize","raster","stringr","velox","assertthat"),
                  .export = c("biodiversity","globalgrid","globalgrid_fraction","output_path",
                              "alignRasters","writeGeoTiff","vx","use_velox")
                  ) %dopar% {
  
        sub <- subset(biodiversity,fname == species)
        stopifnot(nrow(sub)==1) # There should be only one multipolygon in here 
        fname = paste0(sub$GROUP,"_",species,".tif")
        if(file.exists(paste0(output_path,"/",fname))) { return() }
          
        myLog("Rasterizing at ",str_to_title(species))
        
        if(use_velox){
          # Use the velox raster package
          sub$dummy <- 1
          vx2 <- vx$copy()
          ## Rasterize, set background to -1
          vx2$rasterize(spdf=sub, field="dummy", background=0,small = TRUE)
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
        removeTmpFiles(.5)
}
saveRDS(o, paste0('/media/martin/data/features/GARDreptilia',"_",resolution,".rds") )

stopCluster(cl)
stop("DONE!")
stopifnot( length(list.files(output_path)) == length(unique(biodiversity$fname)) )


