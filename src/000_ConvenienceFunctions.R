# Broad functions that most scripts depend on
myLog <- function(...) {
  cat(paste0("[Processing] ", Sys.time(), " | ", ..., "\n"))
}
# Inverse not in function
`%notin%` = function(a, b){!(a %in% b)}

# Normalize a value to a range of 0-1
normalize <- function(x)
{
  if (length(x) == 1) 
    return(1)
  y <- suppressWarnings((x - min(x, na.rm = TRUE))/(max(x, 
                                                        na.rm = TRUE) - min(x, na.rm = TRUE)))
  if (any(is.nan(y))) 
    return(x)
  else return(y)
}

#' Show difference between two vectors
#' 
#' returns a list of the elements of x that are not in y 
#' and the elements of y that are not in x (not the same thing...)
#' @param x vector x
#' @param y vector y
#' @keywords difference
#' @export
#' @usage setdiff2(x,y)
#' @author Martin Jung
setdiff2 <- function(x,y) {
  Xdiff = setdiff(x,y)
  Ydiff = setdiff(y,x)
  list(X_not_in_Y=Xdiff, Y_not_in_X=Ydiff)
}

#### Loading outputs ####
# Load in function
makeStack <- function(files, pol = NULL){
  if(is.null(pol)){
    return( raster::stack(files) )
  } else {
    assertthat::assert_that("sf" %in% class(pol) )
    # Run in parallel
    cl <- registerDoParallel(cores = 7)
    res <- foreach(f = iter(files),.inorder = TRUE,.final = stack) %dopar% {
      vx <- velox(f)
      vx$crop(pol)
      return( vx$as.RasterLayer() )
    }
    n <- str_split(basename(files),"_",3,simplify = T)
    names(res) <- paste(n[,2],str_split(n[,3],"\\.",simplify = T)[,1],sep = "_")
    # Return the result
    return(res)
  }
}

#### Saving outputs ####
# Save output as compressed geotiff file
writeGeoTiff <- function(file, fname,dt = "INT2S"){
  require(assertthat)
  if(!has_extension(fname,"tif")) fname <- paste0(fname,".tif")
  writeRaster(file,fname,
              format='GTiff',
              datatype = dt,
              NAflag = -9999,
              options=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
              overwrite= TRUE )
}

#### Spatial unit conversions ####
latlong2UTMzone <- function(lon,lat){
  # Normal calculation
  ZoneNumber = floor((lon + 180)/6) + 1
  
  # Special case for higher longitude levels
  if( lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0 ){
    ZoneNumber = 32
  }
  # Special cases for svalbard
  if( lat >= 72.0 && lat < 84.0 ) {
    if  ( lon >= 0.0  && lon <  9.0 ) ZoneNumber = 31 
    else if( lon >= 9.0  && lon < 21.0 ) ZoneNumber = 33
    else if( lon >= 21.0 && lon < 33.0 ) ZoneNumber = 35
    else if( lon >= 33.0 && lon < 42.0 ) ZoneNumber = 37
  }
  
  #Return the result
  return(ZoneNumber)
}

DegreeCellAreaKM <- function(lat, height, width) {
  # Returns the area in km squared of a grid cell in degrees of arc
  # lat - the latitudinal centre of the cell
  # height, width - the size of the grid cell in degrees
  
  radians <- function(theta) theta*pi/180.0
  
  # Convert the latitude into radians
  lat.rad <- radians(lat)
  
  # The equatorial and polar radii of the Earth in km
  eq.radius <-  6378137
  pol.radius <- 6356752.3142
  
  # Calculate cell area
  angular.eccentricity <- acos(radians(pol.radius/eq.radius))
  ecc.sq <- sin(radians(angular.eccentricity))^2
  flattening <- 1-cos(radians(angular.eccentricity))
  temp.val <- (eq.radius*cos(lat.rad))^2+(pol.radius*sin(lat.rad))^2
  m.phi <- ((eq.radius*pol.radius)^2)/(temp.val^1.5)
  n.phi <- (eq.radius^2)/sqrt(temp.val)
  lat.length <- pi/180*m.phi/1000
  long.length <- pi/180*cos(lat.rad)*n.phi/1000
  return (lat.length*height*long.length*width)
}

#### Template creation ####
# Empty raster function based on a template
emptyraster <- function(x, ...) { # add name, filename, 
  
  emptyraster <- raster(nrows=nrow(x), ncols=ncol(x),
                        crs=x@crs, 
                        ext=extent(x), ...)
  
  return(emptyraster)
}

#' Unify extent of a rasterstack
#' From Icarus package
#' @author Martin Jung
#' @param rList a list of raster objects

#' @return returns the unified total extent
#' @export
#' 
max_extent <- function(rlist){
  # given list of rasters
  # returns union of extent
  xmin=min(sapply(rl,FUN=function(x){extent(x)@xmin}))
  xmax=max(sapply(rl,FUN=function(x){extent(x)@xmax}))
  ymin=min(sapply(rl,FUN=function(x){extent(x)@ymin}))
  ymax=max(sapply(rl,FUN=function(x){extent(x)@ymax}))
  
  extent(c(xmin,xmax,ymin,ymax))
}


#### Raster Modifications ####
# Define a function for minmax-normalization
normalize_raster <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  return((x - min) / (max - min))
}

#' Allign raster data by bringing it in the same geometry and extend.
#' 
#' If the data is not in the same projection as the template, the allignment
#' will be computed by reprojection only. If the data has already the same
#' projection, the data set will be croped and aggregated prior to resampling
#' in order to reduce computation time.
#' 
#' @author Martin Jung
#' @author Thomas Nauss
#' @import raster
#' @param data raster layer to be resampled
#' @param template raster or spatial data set from which geometry can be extracted
#' @param method method for resampling ("ngb" or "bilinear")
#' @param func Function to apply for aggregation
#' @param cl Should cluster computation be used (Default=T)
#' @param expand Expanding the output variable
#'
#' @return raster layer containing geometrically alligned data
#' @export

alignRasters <- function(data, template, method = "bilinear",func = mean,cl = T,expand = TRUE){
  lib <- c("raster")
  sapply(lib, function(...) stopifnot(require(..., character.only = T)))
  if(cl) beginCluster(parallel::detectCores()-1)
  if(projection(data) == projection(template)){
    data <- raster::crop(data, template, snap = "near")
    if(class(template) == "RasterLayer"){
      if(data@ncols / template@ncols >= 2){
        factor <- floor(data@ncols/template@ncols)
        data <- raster::aggregate(data, fact = factor, fun = func,na.rm = TRUE, 
                          expand=expand)
      }
      data <- raster::resample(data, template, method = method)
    }
  } else {
    data <- projectRaster(data,crs = template, method = method)
  }
  if(cl) endCluster()
  return(data)
}

# Splits a raster into several subtiles to be processed individually
#' @param infile The raster file
#' @param outpath the output path
#' @param parts The number of individual parts used for processing
#' @author Martin Jung
split_raster <- function(infile,outpath,parts=3) {
  ## parts = division applied to each side of raster, i.e. parts = 2 gives 4 tiles, 3 gives 9, etc.
  lib <- c("gdalUtils","parallel","reshape2")
  sapply(lib, function(...) stopifnot(require(..., character.only = T)))
  # Check if infile exists
  stopifnot(file.exists(infile))
  
  filename <- strsplit(basename(infile),".tif")[[1]]
  # Get the dimensions of the file  
  dims <- as.numeric(
    strsplit(gsub('Size is|\\s+', '', grep('Size is', gdalinfo(infile), value=TRUE)), 
             ',')[[1]]
  )
  
  # Generate window frames
  xy <- list()
  # t is nr. of iterations per side
  t <- parts - 1
  for (i in 0:t) {
    for (j in 0:t) {
      # [-srcwin xoff yoff xsize ysize] src_dataset dst_dataset
      srcwin <- paste(i * dims[1]/parts, j * dims[2]/parts, dims[1]/parts, dims[2]/parts)
      xy[[paste0(i,"_",j)]] <- data.frame(infile = infile,srcwin = srcwin, file = paste0(outpath,"/",filename,"_",i,"_",j,".tif"))
    }
  }
  df <- melt(xy)
  
  # Then process per src_win
  cat("Start splitting: ",filename)
  
  # Create a function to split the raster using gdalUtils::gdal_translate
  split <- function(input, outfile, srcwin) {
    gdal_translate(input, outfile, srcwin=srcwin)
  }
  
  # Make a copy for export
  df_org <- df
  # Kick out files already existing
  df <- df[which(!file.exists(as.character(df$file))),]
  
  # Make cluster
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c('split', 'df')) 
  clusterEvalQ(cl,library(gdalUtils))
  
  system.time({
    parLapply(cl, seq_len(nrow(df)), function(i) {
      split(df$infile[i], df$file[i], df$srcwin[i])  
    })
  })
  stopCluster(cl)
  cat("\n")
  cat("Done")
  return(df_org)
}

#' Expand an extent by a certain number
#' 
#' @author Martin Jung
#' @param e an extent object
#' @param f value to increase the extent (Default = 0.1)
#' @return returns the unified total extent
#' @export

extent_expand <- function(e,f=0.1){
  
  xi <- (e@xmax-e@xmin)*(f/2)
  yi <- (e@ymax-e@ymin)*(f/2)
  
  xmin <- e@xmin-xi
  xmax <- e@xmax+xi
  ymin <- e@ymin-yi
  ymax <- e@ymax+yi
  
  return(extent(c(xmin,xmax,ymin,ymax)))  
}

# Transform to mollweide function
projectMollWeide <- function(image, output = "test.tif",s.crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"){
  out <- gdalwarp(srcfile = image@file@name,
                  dstfile = output,
                  s_srs = proj4string(image),
                  t_srs = s.crs,
                  r = "near",
                  co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
                  output_Raster=TRUE, multi = TRUE,
                  overwrite=TRUE,verbose=TRUE)
  return( out )
}

#### Priority map aggregations ####

#' Aggregate a given raster by rank
#' @param ras Calculate the shortfall for each species
#' @param ties.method Default:random
#' @param n The number of breaks
#' @param plot If true, plot the proportions of each cell
#' @import assertthat
#' @import raster
#' @export
#' @author Martin Jung
# New ranked aggregation
raster_to_ranks <- function(ras,ties.method = 'random', n = 101,plot = TRUE){
  assert_that(
    cellStats(ras,'max')>0
  )
  # Make a copy of an empty raster
  mcp_copy <- emptyraster(ras)
  # Convert to data.frame
  df <- as.data.frame(ras)
  # Assign ranks to all values
  df$neworder[!is.na(df[,1])] <- rank(df[!is.na(df[,1]),1],ties.method = ties.method)
  df$newcut <- ggplot2::cut_interval(df$neworder, n = n)
  mcp_copy[] <- df$newcut
  # Subtract 1 from the number of bins to that ranges goes from 0 to 100
  mcp_copy <- mcp_copy - 1
  if(plot){
    # Check that equal number of cells are included in the ranking
    x = cut(mcp_copy[],seq(0,100,10),include.lowest=F)
    #x = cut_number(mcp_ilp_hier[],n = 10)
    barplot(table(x),main = 'Binned by 10%',ylab = 'Number of grid cells',xlab = 'Bins')
  }
  return(mcp_copy)
  # These should be all identical
  #round((table(df$newcut) / length(unique(df$neworder))),2)
}



#### PrioritizeR helper functions ####

# Calculate species-specific relative targets 
# As taken from Mogg et al.
# Resolution being a number in kilometer units
calc_targets <- function(rij_data,resolution) {
  # Resolution is in km!
  rij_data %>% 
    dplyr::group_by(species) %>% # Per species
    dplyr::summarise(
      size = (n_distinct(pu) * (resolution^2)), # Approximate grid cells area occupied in sqkm
      absolute_target = min( c( max( c(2200 ), 0.8 * ((n_distinct(pu) * (resolution^2))  ) ), 10^6 ) )
    ) %>% 
    # Add relative target as well
    dplyr::mutate( relative_target = absolute_target / size ) %>% 
    ungroup()
}

# Alternative implementation using the directly measured area
calc_targets2 <- function(rij_data, multiplier = 1) {
  assert_that(has_name(rij_data,'area'),has_name(rij_data,'species'))
  if(multiplier == 1) { warning('Normal multiplier used. Correct?')}
  # Resolution is in km!
  rij_data %>% 
    mutate(amount = amount * area) %>% 
    dplyr::group_by(species) %>% # Per species
    dplyr::summarise(
      size = sum(amount), # Approximate grid cells area occupied in sqkm
      absolute_target = min( c( max( c(2200 * multiplier ), 0.8 * size  ) ), 10^6 * multiplier )
    ) %>% 
    ungroup() %>% 
    # Add relative target as well
    dplyr::mutate( relative_target = ifelse(  (absolute_target / size)>=1 ,1,(absolute_target / size) ) ) %>% 
    ungroup()
}

# Create a output raster from a prioritizR solution
createOutput <- function(pu_id_raster, r, scol = 'solution_1' ){
  require(assertthat)
  assert_that(has_name(r,scol))
  # Solved object columns: (id, cost, solution_1)
  # Create output raster
  temp <- as.data.frame(pu_id_raster)
  # Replace
  temp$solution <- NA
  temp$solution[ which(temp[,1] %in% r[which(r[,scol] > 0),1] ) ] <- r[which(r[,scol] > 0),scol]
  assert_that(!all(is.na(temp$solution)),msg = 'No solution in output raster')
  
  # Create output template and replace
  solution_raster_data <- emptyraster(pu_id_raster)
  solution_raster_data[] <- temp$solution
  return(solution_raster_data)
}

# Alternative to the feature_representation script
# Necessary as the prioritizeR version creates a M x N matrix
# that is too big for memory (requires ~480GB) for Matrix multiplication
# Here parallizing the process
# #id / cost / protected / solution_1
feature_representation2 <- function(x, solution, cores = 1) {
  if(cores>1){ require(doParallel)}
  require(Matrix);require(progress)
  # Check if any solutions are below 0 (???)
  # TODO: Check with Jeff what this is about? Probably optimization errors/ infesability?
  if(min(unlist(solution), na.rm = TRUE) < 0){
    solution[which(unlist(solution) < 0),] <- 0 # Remove from the solution pool
  }
  # assert valid arguments
  assertthat::assert_that(
    is.data.frame(solution),
    number_of_zones(x) == ncol(solution),
    number_of_total_units(x) == nrow(solution),
    is.data.frame(x$data$cost),
    is.numeric(unlist(solution)),
    min(unlist(solution), na.rm = TRUE) >= 0
#    max(solution, na.rm = TRUE) <= 1 # This was causing the calc to fail?
    )
  # subset planning units with finite cost values
  pos <- x$planning_unit_indices()
  solution <- as.matrix(solution)
  if (any(solution[setdiff(seq_len(nrow(solution)), pos), ,
                   drop = FALSE] > 0))
    stop("planning units with NA cost data have non-zero allocations in the ",
         "argument to solution")
  solution <- solution[pos, , drop = FALSE]
  # calculate amount of each feature in each planning unit
  total <- x$feature_abundances_in_total_units()
  
  if(cores > 1){ registerDoParallel(cores) }
  i=1 #FIXME: Currently only single zone being supported
  # # Calculate blocks
  ll <- seq(1,nrow(x$data$rij_matrix[[i]]) )
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) # Chunking script
  blocks <- chunk(ll, round(length(ll)/ 500) )
  
  if(cores > 1) {
    held <- foreach(bl = iter(blocks),
                    .export = c('x','i','solution'),
                    .packages = c('Matrix'),
                    .inorder = TRUE) %dopar% {
                      return(
                        rowSums( x$data$rij_matrix[[i]][bl,] *  Matrix::Matrix(data = solution[, i],
                                                                             nrow = nrow(x$data$rij_matrix[[i]][bl,]), ncol =  nrow(solution),
                                                                             byrow = TRUE, sparse = TRUE
                          ) 
                        )
                      )
                    }
    # Combine the results
    held <- do.call(c,held)  
    stopImplicitCluster()
  } else {
    held <- vector()
    pb <- progress_bar$new(total = length(blocks))
    for(bl in blocks){
      held <- c(held,
                rowSums( x$data$rij_matrix[[i]][bl,] *  Matrix::Matrix(data = solution[, i],
                                                                       nrow = nrow(x$data$rij_matrix[[i]][bl,]), ncol =  nrow(solution),
                                                                       byrow = TRUE, sparse = TRUE
                  ) 
                )
      )
      pb$tick()
      gc()
    }
    pb$terminate()
    
  }
  stopifnot(
    assertthat::assert_that( length(total) == length(held) )
  )
  out <- tibble::tibble(feature = rep(x$feature_names(), x$number_of_zones()),
                        absolute_held = c(held),
                        relative_held = c(held / total))
  if (x$number_of_zones() > 1) {
    out$zone <- rep(x$zone_names(), each = x$number_of_features())
    out <- out[, c(1, 4, 2, 3), drop = FALSE]
  }
  out
}

## Representation by PU script #
# Script to calculate mean relative representation by PU-ID
# Output is meant to create a map showing the irreplaceability of each pu
feature_pu_representation <- function(x,solution,out,pu_id_raster){
  require(doParallel);require(Matrix);require(progress)
  # assert valid arguments
  assertthat::assert_that(
    is.data.frame(solution),
    number_of_zones(x) == ncol(solution),
    number_of_total_units(x) == nrow(solution),
    is.data.frame(x$data$cost),
    is.numeric(unlist(solution)),
    min(unlist(solution), na.rm = TRUE) >= 0,
    max(unlist(solution), na.rm = TRUE) <= 1,
    "absolute_held" %in% names(out)
    )
  # subset planning units with finite cost values
  pos <- x$planning_unit_indices()
  solution <- as.matrix(solution)
  if (any(solution[setdiff(seq_len(nrow(solution)), pos), ,
                   drop = FALSE] > 0))
    stop("planning units with NA cost data have non-zero allocations in the ",
         "argument to solution")
  solution <- solution[pos, , drop = FALSE]
  # Assign NA data to non-represented features
  out[which(out[,'absolute_held']==0),'absolute_held'] <- NA
  
  i=1
  ll <- seq(1,nrow(solution) )
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) # Chunking script
  blocks <- chunk(ll, round(length(ll)/ 500) )
 
  held <- vector()
  pb <- progress_bar$new(total = length(blocks))
  for(bl in blocks){
    # Multiply by total feature abundance divisor
    held <- c( held,
              colSums((x$data$rij_matrix[[i]][,bl] * (1/out[,'absolute_held'])) * Matrix::Matrix(data = solution[bl, i],
                                                                                      nrow = nrow(x$data$rij_matrix[[i]]), ncol =  length(solution[bl,]),
                                                                                      byrow = TRUE, sparse = TRUE
              ),na.rm = TRUE)
             )
    pb$tick()
  }
  pb$terminate()
  # Securitz checks
  stopifnot(
    assertthat::assert_that( 
      nrow(solution) == length(held)
      )
  )
  # Create output raster
  temp <- as.data.frame(pu_id_raster)
  # Replace
  temp$solution <- NA
  temp[match(pos,temp[,1]),"solution"] <- held
  temp$solution[which(temp$solution==0)] <- NA # Sum of Zero is not relevant
  assert_that(!all(is.na(temp$solution)),msg = 'No solution in output raster')
  
  # Create output template and replace
  solution_raster_data <- emptyraster(pu_id_raster)
  solution_raster_data[] <- temp$solution
  solution_raster_data <- normalize_raster(solution_raster_data) # Normalize
  names(solution_raster_data) <- "complementarity"
  #writeGeoTiff(solution_raster_data, paste0(results_path,'/','maxfeat_speciestargets_complementarity_',target_resolution,'_',round(b,2)*100,'perc.tif'),'FLT4S')
  return(solution_raster_data)
}

#' Calculate shortfall for a given representation file
#' Possible either by species or overall 
#' @param repr_file Calculate the shortfall for each species
#' @param by_species Calculate shortfall per species or how many targets have been reached
#' @export
#' @author Martin Jung

calc_shortfall <- function(repr_file, by_species = TRUE){
  require(assertthat):require(tidyverse)
  assert_that(
    is.data.frame(repr_file),is.logical(by_species),
    has_name(repr_file, "absolute_held"),has_name(repr_file, "absolute_target")
  )
  if(by_species){
    # Calculate the shortfall for each feature
    rr <- repr_file %>% mutate(
      shortfall = (pmax(absolute_target - absolute_held, 0))/ absolute_target,
      target_reached = absolute_held >= absolute_target
    ) %>%
      dplyr::select(feature,shortfall,target_reached)
  } else {
    # Simply calculate overall how many targets have been reached
    rr <- repr_file %>% 
      dplyr::summarise(target_reached = length(which( absolute_held >= absolute_target )) / n() ) %>% 
      dplyr::mutate(feature = 'OVERALL') %>% dplyr::select(feature, target_reached)
  }
  return(rr)
}