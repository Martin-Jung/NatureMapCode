# This is the central figure script to create final output figures 
# for all global and sub analyses
# Selection of parameters at the start specifies which files are to be loaded

## Paths
results_path <- 'results' # Results path
figure_path <- 'figures' # Figure output path
output_path <- 'data' # Data path
temp_path <- tempdir()
if(!dir.exists(figure_path)) dir.create(figure_path)

target_resolution <- "10km";target_resolution_number <- as.numeric(gsub("\\D","",target_resolution))
projection <- "mollweide"

rij_path <- paste0(output_path,"/table_puspecies",'_',target_resolution,'.fst') # The final output table
feature_path <- paste0(output_path,"/speciesID_table_",target_resolution,".fst") # The table for feature preparation
biome_path <- paste0(output_path,"/pu_biome_",target_resolution,".fst")
path_data <- '/media/martin/data/'
ecoregion_path <- paste0(path_data,"raw/Ecoregions2017/Ecoregions2017.shp")
gadm_path <- paste0(path_data,"raw/gadm36_shp")
carbon_path <- paste0("features_esh/carbon_agbc")
water_path <- paste0("features_esh/water")

repr_id <- TRUE # Where representative sets used?
repr_run_path <- paste0('REPrun_') # The path for  the representative sets

# Output data formats
target_range <- c("","_esh")[2] # Range or ESH?
curr = "15perc";curr_number <- 15 # What is the current data if existing
pa_fname <- c("",'withPA')[2] # Solution with PA's excluded ? 
pa_carbon <- c("",'_carbon_')[1] # Carbon be included as linear feature ?
pa_water <- c("","_water_")[1] # Water included
split_id <- c("","_biome.id_")[2] # Was a biome stratification used?
pa_carbmult <- c("","_carbweight2_")[1] # Differing carbon weights
pa_watermult <- c("","_waterweight2_")[1] # Different water weight
pa_plants <- c("_waPlants","")[2] # Plants included?
pa_phylo <- c('','_phylo-rank_weight-comparison_','_phylo-ed-comparison_',
              '_phylo-edge-comparison_','_phylo-noweight-comparison_')[1] # Phylogenetic weights included

# ------------------------ #
# Code start #
library(ggplot2)
library(colorspace);library(scico);library(wesanderson)
library(extrafont)
library(ggrepel)
library(ggridges)
library(scales)
library(colorspace)
library(raster)
library(gdalUtils)
library(sf)
library(ggthemes)
library(tidyverse)
library(fst)
library(pbapply)
library(RStoolbox)
library(tmap);library(tmaptools)
library(assertthat)
source('src/000_ConvenienceFunctions.R')
stopifnot(
  assert_that(
    file.exists(ecoregion_path),
    dir.exists(output_path),
    dir.exists(results_path),
    file.exists(feature_path),
    exists("makeStack")
  )
)
# Land area
land <- raster(paste0("data/globalgrid_",projection,"_",target_resolution,".tif"))
land <- setMinMax(land)
land.mask <- land
land.mask[land.mask>0] <- 1;land.mask[land.mask!=1] <- NA
# ------------------------ #
#### [Figure] Global ranked asset mapping ####
# Organized part to rank all computed assets globally
# Get budgets

if(repr_id){
  dd <- list.dirs(paste0(results_path,"/",target_resolution,target_range,"/"))
  dd <- dd[grep(repr_run_path,dd)]
  # Get list of representative files
  ll <- list.files(path = dd,
                   pattern = paste0(
                     'minshort_speciestargets',
                     split_id,
                     ifelse(pa_fname=="","",pa_fname),
                     pa_carbon,
                     pa_water,
                     pa_carbmult,
                     pa_watermult,
                     pa_plants,
                     pa_phylo,
                     target_range,
                     target_resolution
                   ),
                   full.names = TRUE
  )
  number_representativesets <- length(dd)
  
  # Get list of files with shortfall
  ll.short <- list.files(path = dd,
                   pattern = paste0(
                     'minshort_speciestargets',
                     split_id,
                     ifelse(pa_fname=="","",pa_fname),
                     pa_carbon,
                     pa_water,
                     pa_carbmult,
                     pa_watermult,
                     pa_plants,
                     pa_phylo,
                     target_range,
                     target_resolution
                   ),
                   full.names = TRUE
  )
  ll.short <- ll.short[has_extension(ll.short,'fst')]

} else {
  ll <- list.files(path = paste0(results_path,"/",target_resolution,target_range,"/"),
                   pattern = paste0(
                     'minshort_speciestargets',
                     split_id,
                     ifelse(pa_fname=="","",pa_fname),
                     pa_carbon,
                     pa_water,
                     pa_carbmult,
                     pa_watermult,
                     pa_plants,
                     pa_phylo,
                     target_range,
                     target_resolution
                   ),
                   full.names = TRUE
  )
}
stopifnot(assert_that(length(ll)>0))

# Get tif files 
ll.tif <- ll[has_extension(ll,"tif")]
(ll.tif)

# Output filename
out_dir <- paste0(results_path,"/",target_resolution,target_range,"_ranked")
if(!dir.exists(out_dir)) dir.create(out_dir)
(fname <- paste0(out_dir,"/minshort_speciestargets",
                 split_id,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,
                 ifelse(repr_id,paste0('_repruns',number_representativesets),""),
                 "_ranked.tif"))
(fname.cv <- paste0(out_dir,"/minshort_speciestargets",
                 split_id,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,pa_phylo,target_range,target_resolution,
                 ifelse(repr_id,paste0('_repruns',number_representativesets),""),
                 "_cv.tif"))

# Biodiversity value(s)
gc(); raster::removeTmpFiles(.2)
# Load individual solutions (raster) from the representative sets
ras2 <- makeStack(ll.tif)
# Reclassifies all PUs in the stack to have 0 instead of NA so that 
# averaging the selected proportion takes account of non-selected cells
ras2 <- reclassify(ras2, cbind(NA, NA, 0), right = FALSE)

# Average stack with all proportions
mcp_ilp_hier <- mean(ras2, na.rm = T) # Mean! 
# Mask with global land area mask
mcp_ilp_hier <- raster::mask(mcp_ilp_hier, land.mask)

# Multiply with 100 (precision of input layers) for the ranking
mcp_ilp_hier_norm <- mcp_ilp_hier * 100 
# Remove 0 fractional estimates so that equal-width bins are constructed correctly
mcp_ilp_hier_norm[mcp_ilp_hier_norm==0] <- NA
# Convert to ranks and then bin (equal-width) averages across PUs. Ocassional ties are resolved randomly
# The resulting ranked map is then inverted so that priority ranks go from 1 to 100
mcp_ilp_hier_norm <- abs( raster_to_ranks(mcp_ilp_hier_norm,
                                          globalgrid = globalgrid,
                                          method = 'area',n = 101,plot = T) # - 100 
)
# land30 <- mcp_ilp_hier_norm
# land30[land30 > 30] <- NA; land30[land30 > 0] <- 1; (cellStats(land30 * globalgrid,'sum') - (cellStats(globalgrid, 'sum') * .3) )/(cellStats(globalgrid, 'sum') * .3)
# Add least important areas as 100 so that proportions match again. This is necessary for prioritization variants where not all PUs have value in reaching the targets
mcp_ilp_hier_norm[is.na(mcp_ilp_hier_norm)] <- 100 
# Mask again with global land area mask 
mcp_ilp_hier_norm <-  raster::mask(mcp_ilp_hier_norm, land.mask)
assertthat::assert_that(cellStats(mcp_ilp_hier_norm,min)==1, cellStats(mcp_ilp_hier_norm,max)==100)
NAvalue(mcp_ilp_hier_norm) <- -9999
proj4string(mcp_ilp_hier_norm) <- proj4string(globalgrid)

# Test plot
plot(mcp_ilp_hier_norm,col = scico::scico(10,palette = 'roma'))

# Save outputs #
writeGeoTiff(mcp_ilp_hier_norm, fname = fname) # Save outputs and clear
# --- #

# For variation across representative sets, use the coefficient of variation
co.var <- function(x, na.rm = TRUE) ( sd(x,na.rm = TRUE) / mean(x,na.rm = TRUE) )
ras3 <- abs((ras2*100)-100)
mcp_covar <- calc(ras3 ,co.var)
# Mask with global land area mask
mcp_covar[is.na(mcp_covar)] <- 0
mcp_covar <-  raster::mask(mcp_covar,land)
NAvalue(mcp_covar) <- -9999
plot(mcp_covar, col = scico(10,palette = 'tokyo',direction = -1))
# Save outputs #
writeGeoTiff(mcp_covar, fname = fname.cv, dt='FLT4S') # Save outputs and clear  

# ----------------------------------------------------------- #
#### [Figure] Nice plotting of full maps for Figure 1 and SI-Figure 1 ####
# ----------------------------------------------------------- #
ras_mean <- raster(fname)
ras_cv <- raster(fname.cv)

(out_fname <- paste0(figure_path, "/",tools::file_path_sans_ext(basename(fname)), ".png" ))

# Nice plotting of the selected priority map
g <- ggR(ras_mean, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_gradientn(breaks = c(100,75,50,25,1),colours = scico(10,palette = 'roma'),na.value = 'transparent',
                       guide = guide_colourbar(title = "Priority rank",title.theme = element_text(size = 14),
                                               title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                               direction = 'horizontal',ticks = F,reverse = T) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "")
ggsave(filename = out_fname, plot = g,width = 8, height = 8, dpi = 400)

# Plot for the uncertainty only
(out_fname_cv <- paste0(figure_path, "/",tools::file_path_sans_ext(basename(fname.cv)), ".png" ))
g <- ggR(ras_cv, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_viridis_c(option = 'D',direction = -1, na.value = 'transparent',breaks = c(0,2.5,5,7.5,10),
                       guide = guide_colourbar(title = "Rank uncertainty (%)",title.theme = element_text(size = 14),
                                               title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                               direction = 'horizontal',ticks = F,reverse = T) ) +
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = out_fname_cv, plot = g,width = 8, height = 8, dpi = 400)

# inverted and as SVG
# Insert radial legend in here
library(ggradar)
stopifnot( assert_that(length(ll.short)>0))
perc10 <- bind_rows(
  pblapply(ll.short[grep('10perc',ll.short)], function(x) {calc_shortfall(read_fst(x)) %>% 
      dplyr::mutate(run = basename(dirname(x)),budget = "10") } ),
  pblapply(ll.short[grep('30perc',ll.short)], function(x) {calc_shortfall(read_fst(x)) %>% 
      dplyr::mutate(run = basename(dirname(x)),budget = "30") } )
)

perc10 <- perc10[-which(is.na(perc10$shortfall)),]

df <- bind_rows(
  perc10 %>% dplyr::filter(feature %notin% c('AGB_carbon','Water')) %>% 
    dplyr::group_by(feature,budget) %>% 
    dplyr::summarise(prop = mean(shortfall, na.rm = TRUE) ) %>% ungroup() %>% 
    dplyr::group_by(budget) %>% 
    dplyr::summarise(prop = sum(prop==0, na.rm = TRUE) / n() ) %>% 
    dplyr::mutate(biocarbongroup = 'Biodiversity'),
  perc10 %>% dplyr::filter(feature %in% c('AGB_carbon','Water')) %>% 
    dplyr::group_by(budget,feature) %>% 
    dplyr::summarise(prop = mean(shortfall) ) %>% 
    dplyr::mutate(feature = factor(feature, levels = c('AGB_carbon','Water'), labels = c('Carbon','Water') )) %>% 
    dplyr::rename(biocarbongroup = feature)
) %>% 
  tidyr::spread(biocarbongroup, prop) %>% dplyr::mutate(group = 1) %>% dplyr::select(group,budget,Biodiversity:Water) %>% 
  mutate(budget = as.numeric(budget))

df %>% dplyr::select(budget,Biodiversity,Carbon,Water) %>% 
  dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% 
  dplyr::mutate(Biodiversity = abs(Biodiversity - 100))

# Colours
cols = c( scico(10,palette = 'roma')[c(1,3)] )

gr <- df %>% dplyr::select(budget,Biodiversity,Carbon,Water) %>% 
  dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% 
  dplyr::mutate(Biodiversity = abs(Biodiversity - 100)) %>% 
  ggradar(.,
          grid.min = 0,grid.mid = 50, grid.max = 100,
          axis.label.size = 6, grid.label.size = 7, gridline.label.offset = 0,legend.text.size = 10,
          base.size = 20, group.point.size = 4,group.line.width = 1.05
          #          axis.labels = c('Biodiversity','Carbon','Water'), plot.title = ''
  ) +
  theme_void(base_size = 8) + 
  theme(legend.title = element_text(size = 12, angle = 90),legend.text = element_text(size = 12),
        legend.position = 'left') + guides(colour = guide_legend(title = "Total land area (%)",nrow = 1)) +
  scale_colour_manual(values = cols) #+ guides(colour = 'none') 
gr
(out_fname <- paste0(figure_path, "/",tools::file_path_sans_ext(basename(fname)), "_spiderplot.png" ))
(out_fname <- paste0(figure_path, "/",tools::file_path_sans_ext(basename(fname)), "_spiderplot.svg" ))
ggsave(filename = out_fname,plot = gr, width = 5, height = 5, dpi = 400,bg = 'transparent')

# -------------------------------- #
#### [Figure] Effect of adding plant species - Shift of global ranks ####
# First create the global map
number_representativesets = 10
out_dir <- paste0(results_path,"/",target_resolution,target_range,"_ranked")
ras_withplants <- raster( 'results/10km_esh_ranked/minshort_speciestargets_biome.id__esh10km_repruns10_ranked.tif')
ras_explants <- raster( 'results/10km_esh_ranked/minshort_speciestargets_biome.id__waPlants_esh10km_ranked.tif' )

# Figure output
(
  out_fname <- paste0(figure_path, "/minshort_speciestargets_waPlants",target_range,target_resolution,
                    ifelse(repr_id,paste0('_repruns',number_representativesets),""),
                    "_ranked_difference.png")
)
# Positive values mean higher priorities
# Negative: Removing plants increased priorities
# Positive: Removing plants decreased priorities
ras_viz <- (ras_explants - ras_withplants) 
writeGeoTiff(ras_viz, fname = paste0(out_dir,paste0("/waplants_",
                                                    ifelse(repr_id,paste0('_repruns',number_representativesets),""),
                                                    "_rank_difference.tif")),'INT2S') # Save outputs and clear 

cols <- scico(n = 7,palette = 'vikO',direction = -1) # cork or broc or lisbon or tofino
cols[4] <- "#F2F2F2"

# The plot
g <- ggR(ras_viz, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_gradientn(colours = cols,na.value = 'transparent',
                       guide = guide_colourbar(title = "Change in\npriority rank",title.theme = element_text(size = 14),
                              title.position = 'bottom',title.hjust = 0,label.position = 'bottom',barwidth = unit(1.5,'in'),
                              direction = 'horizontal',ticks = F,reverse = F) ) +
  theme(legend.position = c(.1,.1),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = out_fname, plot = g,width = 10, height = 14, dpi = 400)

# Also make a histogram of t  he differences
vals <- data.frame(values = ras_viz[]) %>% drop_na()
gg <- ggplot(vals,aes(x = values)) +
        theme_classic(base_size = 20,base_family = 'Arial') +
        geom_histogram(bins = 50,fill = NA,colour = 'black') +
  scale_y_continuous(expand = c(0,0),limits = c(0,220000),position = 'left') +
  labs(x = "Change in priority rank", y = 'Frequency (10km² grid cells)')
gg
ggsave(filename = "figures/minshort_speciestargets_waPlants_esh10km_repruns10_ranked_difference_histogram.png",
       plot = gg,width = 6, height = 5, dpi = 400,bg = 'transparent')

# --- #
# Next calculate species representation curves with and without plants
target_resolution <- '10km'
pa_fname = ''
split_id = '_biome.id_'
dd <- list.dirs(paste0(results_path,"/",target_resolution,target_range,"/"))
if(repr_id) { dd <- dd[grep(repr_run_path,dd)] } else { dd <- dd[1]}
list_withplants <- list.files(path = dd,
                                    pattern = paste0(
                                      'minshort_speciestargets',
                                      split_id,
                                      "",
                                      "", #pa_carbon,
                                      "", #pa_water,
                                      pa_carbmult,
                                      pa_watermult,
                                      "",
                                      pa_phylo,
                                      target_range,
                                      target_resolution
                                    ),
                                    full.names = TRUE
)
list_withplants <- list_withplants[has_extension(list_withplants,"fst")]
list_withoutplants <- list.files(path = list.dirs(paste0(results_path,"/",target_resolution,target_range,"/"))[1],
                                       pattern = paste0(
                                         'minshort_speciestargets',
                                         split_id,
                                         "",
                                         "", #pa_carbon,
                                         "", #pa_water,
                                         pa_carbmult,
                                         pa_watermult,
                                         "_waPlants", #pa_plants,
                                         pa_phylo,
                                         target_range,
                                         target_resolution
                                       ),
                                       full.names = TRUE
)
list_withoutplants <- list_withoutplants[has_extension(list_withoutplants,"fst")]

stopifnot(assert_that(length(list_withplants) > 0, length(list_withoutplants) > 0 ))
# All representation budgets
rr <- bind_rows(
    do.call(rbind, pblapply(list_withplants, function(x,nr = 6) { #4
      read.fst(x) %>% calc_shortfall %>%
        dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
    }) ) %>% dplyr::mutate(type = "Vertebrates and plants"),
    do.call(rbind, pblapply(list_withoutplants, function(x,nr = 7) { #5
      read.fst(x) %>% calc_shortfall %>%
        dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
    })
  ) %>% dplyr::mutate(type = "Vertebrates only")
)
rr$budget <- gsub('\\D','',rr$budget)
rr$budget <- factor(rr$budget, levels = c(if(pa_fname!="") {15}, seq(ifelse(pa_fname=="",10,20),100,10)) )
rr <- rr[-which(is.na(rr$shortfall)),]
# Assure that only vertebrate species are included and identical in both
vertebrate_species <- rr %>% filter(type == "Vertebrates only") %>% select(feature) %>% distinct()
rr <- rr %>% filter(feature %in% vertebrate_species$feature)

# Summarize
rr_full <- rr %>% dplyr::group_by(type, run,budget) %>% 
  dplyr::summarise(target_reached = sum(target_reached,na.rm = TRUE) / n() ) %>% ungroup() %>% 
  dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
  # Average + sd per run
  dplyr::group_by(type,budget) %>% 
  dplyr::summarise(target_reached_avg = mean(target_reached*100,na.rm=TRUE),
                   target_reached_ymin = mean(target_reached*100,na.rm=TRUE) - sd(target_reached*100,na.rm=TRUE),
                   target_reached_ymax = mean(target_reached*100,na.rm=TRUE) + sd(target_reached*100,na.rm=TRUE))

# Plot for plants only
gp <- rr_full %>% 
  # ----- #
  #dplyr::filter(budget <= 60) %>% 
  ggplot(., aes(x = budget, y = target_reached_avg, group = type, linetype = type)) + 
  theme_few(base_size = 20) +
  #geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
  geom_line(size = 2,alpha = .65) +
  #geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax), fill = 'grey80', alpha = .1) +
  scale_linetype_manual(values = c('solid','twodash'),guide = guide_legend("")) + 
  theme(legend.position = c(.65,.15),legend.key.width = unit(.5,'in'),legend.key = element_rect(fill = 'transparent') ) +
  scale_colour_manual(values = 'black') + guides(colour = 'none', fill = 'none') +
  scale_x_continuous(breaks = pretty_breaks(5),limits = c(as.numeric(levels(rr$budget)[1]),100)) +
  scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
  ylab("Targets reached (% of all assets)") + xlab("Area managed for conservation\n(% of land area)") +
  theme(plot.background = element_blank())
gp
ggsave(filename = paste0(figure_path, "/minshort_speciestargets_waPlants",target_range,target_resolution,"_repruns10_representation_difference.png"),
       plot = gp, width = 6, height = 6,dpi = 400)
# Stats
rr_full %>% dplyr::filter(budget == '10')
rr_full %>% select(type, budget,target_reached_avg) %>% 
  tidyr::pivot_wider(id_cols = 'budget',names_from = 'type', values_from = 'target_reached_avg') %>% 
  mutate(diff = (`Vertebrates and plants` - `Vertebrates only`))

# Which biomes overall ?
biomes <- raster(paste0( '/media/martin/data/constraints/biomes_',target_resolution,'_',projection,'_modal.tif') )
biomes <- alignRasters(biomes, ras_viz, method = 'ngb', func = raster::modal, cl = TRUE)
stopifnot( assert_that( compareRaster(ras_viz,biomes) ) )
# The biomes colour code for plotting
biomes_code <- sf::read_sf('/media/martin/data/raw/Ecoregions2017/Ecoregions2017.shp') %>% 
  dplyr::select(BIOME_NUM,BIOME_NAME,COLOR_BIO) %>% sf::st_drop_geometry() %>% distinct() %>% 
  dplyr::slice(-15) %>% 
  dplyr::arrange(BIOME_NUM)
biomes_code$BIOME_NAME <- factor(biomes_code$BIOME_NAME,levels = biomes_code$BIOME_NAME) # Convert to factor using current order

# Recode the biome names
# Calculate the mean difference at each 
df <- data.frame(biome = getValues(biomes),
                 diff = getValues(ras_viz)
) %>% drop_na() %>% 
  # Summarise
  dplyr::group_by(biome) %>% 
  dplyr::summarise(diff_mean = mean(diff),
                   diff_ymin = mean(diff) - sd(diff),
                   diff_ymax = mean(diff) + sd(diff)
                   )

# Left_join with calculation
df <- left_join(df, biomes_code, by = c('biome' = 'BIOME_NUM'))

# Make the plot
gb <- ggplot(df, aes(x = fct_reorder(BIOME_NAME,diff_mean), y = diff_mean)) +
  theme_few(base_size = 20) + 
  coord_flip() +
  geom_point(size = 2, colour = 'black') +
  geom_pointrange(aes(ymin = diff_ymin, ymax = diff_ymax), size = 1.5, stroke = 1.5) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  #scale_colour_manual(values = df$COLOR_BIO,guide = guide_legend(title = "",ncol = 3)) +
  scale_y_continuous(breaks = pretty_breaks(5),expand = c(0,0)) +
  scale_x_discrete(position = "top") + 
  theme(axis.text.y = element_text(size = 18)) +
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(legend.position = 'bottom',legend.text = element_text(size = 10)) +
  labs(x = "", y = "Change in priority rank") +
  guides(colour = 'none')
gb
ggsave(filename = paste0(figure_path, "/minshort_speciestargets_waPlants",target_range,target_resolution,"_repruns10_biome_difference.png"),
       plot = gb, width = 12, height = 6,dpi = 400)

# --------------------- #
#### [Figure] Representation analysis per budget and groups of species ####
# Make plots for the representation
# Load all the files of the given solution
myLog("Loading results now with increasing budget")
# Load  
target_resolution <- '10km'
pa_fname <- c("","withPA")[1] # Locked in
pa_carbon <- c("",'_carbon_')[1] # Carbon be included as linear feature ?
pa_water <- c("","_water_")[1] # Water included
split_id <- c("","_biome.id_")[1]
dd <- list.dirs(paste0(results_path,"/",target_resolution,target_range,"/"))
if(repr_id) { dd <- dd[grep(repr_run_path,dd)] } else { dd <- dd[1]}
# Get list of files with shortfall
ll <- list.files(path = dd,
                       pattern = paste0(
                         'minshort_speciestargets',
                         split_id,
                         ifelse(pa_fname=="","",pa_fname),
                         pa_carbon,
                         pa_water,
                         pa_carbmult,
                         pa_watermult,
                         pa_plants,
                         pa_phylo,
                         target_range,
                         target_resolution
                       ),
                       full.names = TRUE
)
ll <- ll[has_extension(ll,'fst')]
stopifnot(length(ll)>0);length(ll)

# All representation budgets
rr <- do.call(rbind, pblapply(ll, function(x) {
      nr = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)
      nr = which(apply(nr, 2, function(y) str_detect(y,'perc')))
      read.fst(x) %>% calc_shortfall %>%
                dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
              })
)
stopifnot( assert_that('100perc' %in% rr$budget) )
# Join in species ids
if(split_id != "_biome.id_"){
  rr <- rr %>% left_join(.,
                         read_fst(feature_path) %>% mutate(feature = paste0(iucn_id_no,'__',binomial) ) %>%
                           dplyr::select(feature,category, data, kingdom, class) %>% 
                           dplyr::bind_rows( tibble()  ), # Add Carbon!
                         by ='feature' )
  stopifnot( assert_that(anyNA(rr$feature)==FALSE) )
} else {
  # In case we use the biome split
  ff <- read_fst(feature_path) %>% mutate(feature = str_to_lower( binomial ) ) %>%
    dplyr::select(iucn_id_no,feature,category, data, kingdom, class) %>% 
    dplyr::bind_rows( tibble() )
  # Get ids from sets and only those species that in the id list
  lf <- list.files('sets','splist',full.names = T)
  repr_id_file <- lapply(lf, function(x) { read_csv(x) %>% dplyr::select(2) } ) %>% 
    bind_rows() %>% distinct()
  ff <- ff %>% filter(iucn_id_no %in% repr_id_file$x) %>% dplyr::select(-iucn_id_no)
  rm(lf, repr_id_file)
  # --- #
  rr <- rr %>% mutate(feature = str_to_lower( str_split(rr$feature,"__",simplify = T)[,1]) ) %>% 
    left_join(.,ff, # Add Carbon!
                by ='feature' )
  stopifnot( assert_that(anyNA(rr$feature)==FALSE) )
}
# Reformat
rr$budget <- gsub('\\D','',rr$budget)
if(pa_fname!=""){ rr$budget[rr$budget %in% c(14,15)] <- "15"} 
rr$budget <- factor(rr$budget, levels = c(if(pa_fname!="") {15}, seq(ifelse(pa_fname=="",10,20),100,10)) )
# Set data for carbon and water
rr[which(rr$feature=="agb_carbon"),c('data')] <- 'Carbon'
rr[which(rr$feature=="water"),c('data')] <- 'Water'
# Make a class for taxonomic grouping
stopifnot( assert_that(!all(is.na(rr$budget)),anyNA(rr$data)==FALSE ) )

# --------------------------------------------------------- #
# Link with taxonomic information
rr <- rr %>% # Assign taxonomic grouping
  mutate(TGrouping = factor(data,
                            levels = c( "IUCN Mammal", "IUCN Reptiles", "IUCN Amphibians", "IUCN Bird", "GARD Shai",
                                        "IUCN Plants", "Kew IUCN Plants", "Kew BGCI Plants", "BIEN Plantae Points", "New Plants - PPM",     
                                        "New Plants - Rangebag","BIEN Plantae PPM", "New Plants - Points",
                                        "Carbon", "Water" 
                                        ),
                            labels = c('Mammals', 'Reptiles', 'Amphibians', 'Birds','Reptiles',
                                       'Plants','Plants','Plants','Plants','Plants','Plants','Plants','Plants',
                                       'Carbon','Water') )
  ) %>% mutate(category = ifelse(is.na(category),"DD",category )) %>%
  mutate(category = fct_collapse(category, DD = c("LR/cd","LR/lc","LR/nt","DD"))) %>% 
  mutate(
    category = factor(category,levels = c("EX","EW","CR","EN","VU","NT","LC","DD"))
  ) %>% mutate(category = droplevels(category)) %>% # Ensure that extinct species are truly removed!
  # Another group for threatened / non-threatened species
  mutate(category.threatened = fct_collapse(category, threatened = c("CR","EN","VU"),
                                            `non-threatened` = c("NT","LC","DD")) 
  ) 

# Remove carbon and water here
rr <- rr %>% dplyr::filter(TGrouping %notin% c('Carbon','Water'))
if(pa_carbon == "_carbon_" && pa_water == "_water_"){
  myLog('Bio,Water and Carbon')
  write_fst(rr,paste0("temporary_results/temporary_biocarbonwater",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
} else 
if(pa_carbon == "_carbon_" && pa_water == "") {
  myLog('BioCarbon only')
  write_fst(rr,paste0("temporary_results/temporary_biocarbon",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
} else
if(pa_water == "_water_" && pa_carbon == ""){
  myLog('BioWater only')
  write_fst(rr,paste0("temporary_results/temporary_biowater",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
} else
{  myLog('Bio only')
  write_fst(rr,paste0("temporary_results/temporary_bioonly",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
}
# ------------------------------ #
rr <- read_fst(paste0("temporary_results/temporary_bioonly",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
rr2 <- read_fst(paste0("temporary_results/temporary_biocarbon",pa_fname,"_",target_resolution,split_id,"_reprsave.fst")) # Get the separate results
rr3 <- read_fst(paste0("temporary_results/temporary_biocarbonwater",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))
rr4 <- read_fst(paste0("temporary_results/temporary_biowater",pa_fname,"_",target_resolution,split_id,"_reprsave.fst"))

# Build overall dataset
rr_full <- bind_rows(
    rr %>% dplyr::group_by(run,feature,budget) %>% 
    dplyr::summarise(target_reached = sum(target_reached,na.rm = TRUE) / n() ) %>% ungroup() %>% 
    dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
    # Average + sd per run
    dplyr::group_by(feature,budget) %>% 
    dplyr::summarise(target_reached = mean(target_reached,na.rm = TRUE)) %>% ungroup() %>% 
    dplyr::group_by(budget) %>% 
    dplyr::summarise(target_reached_avg = mean(target_reached*100,na.rm = TRUE)
#                     target_reached_ymin = mean(target_reached*100,na.rm = TRUE) - sd(target_reached*100,na.rm = TRUE),
#                     target_reached_ymax = mean(target_reached*100,na.rm = TRUE) + sd(target_reached*100,na.rm = TRUE)
                     ) %>% 
    dplyr::mutate(type = 'Biodiversity only'),
    rr2 %>% dplyr::group_by(run,feature,budget) %>% 
      dplyr::summarise(target_reached = sum(target_reached,na.rm = TRUE) / n() ) %>% ungroup() %>% 
      dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
      # Average + sd per run
      dplyr::group_by(feature,budget) %>% 
      dplyr::summarise(target_reached = mean(target_reached,na.rm = TRUE)) %>% ungroup() %>% 
      dplyr::group_by(budget) %>% 
      dplyr::summarise(target_reached_avg = mean(target_reached*100,na.rm = TRUE)
#                       target_reached_ymin = mean(target_reached*100,na.rm = TRUE) - sd(target_reached*100,na.rm = TRUE),
#                       target_reached_ymax = mean(target_reached*100,na.rm = TRUE) + sd(target_reached*100,na.rm = TRUE)
                       ) %>% 
      dplyr::mutate(type = 'Biodiversity and Carbon'),
    rr4 %>% dplyr::group_by(run,feature,budget) %>% 
      dplyr::summarise(target_reached = sum(target_reached,na.rm = TRUE) / n() ) %>% ungroup() %>% 
      dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
      # Average + sd per run
      dplyr::group_by(feature,budget) %>% 
      dplyr::summarise(target_reached = mean(target_reached,na.rm = TRUE)) %>% ungroup() %>% 
      dplyr::group_by(budget) %>% 
      dplyr::summarise(target_reached_avg = mean(target_reached*100,na.rm = TRUE)
#                       target_reached_ymin = mean(target_reached*100,na.rm = TRUE) - sd(target_reached*100,na.rm = TRUE),
#                       target_reached_ymax = mean(target_reached*100,na.rm = TRUE) + sd(target_reached*100,na.rm = TRUE)
                       ) %>% 
      dplyr::mutate(type = 'Biodiversity and Water'),
    rr3 %>% dplyr::group_by(run,feature,budget) %>% 
      dplyr::summarise(target_reached = sum(target_reached,na.rm = TRUE) / n() ) %>% ungroup() %>% 
      dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
      # Average + sd per run
      dplyr::group_by(feature,budget) %>% 
      dplyr::summarise(target_reached = mean(target_reached)) %>% ungroup() %>% 
      dplyr::group_by(budget) %>% 
      dplyr::summarise(target_reached_avg = mean(target_reached*100,na.rm = TRUE)
#                       target_reached_ymin = mean(target_reached*100,na.rm = TRUE) - sd(target_reached*100,na.rm = TRUE),
#                       target_reached_ymax = mean(target_reached*100,na.rm = TRUE) + sd(target_reached*100,na.rm = TRUE)
                       ) %>% 
      dplyr::mutate(type = 'Biodiversity, Carbon and Water')
)
rr_full$type <- factor(rr_full$type, levels = c("Biodiversity only", "Biodiversity and Carbon","Biodiversity and Water", "Biodiversity, Carbon and Water"))

rr_full_forexport <- rr_full

# -------------------------- #
# Some summary stats - Shortfall at given budgets or taxonimic groups

# Overall shortfall for all assets
stopifnot(assert_that('Carbon' %in% unique(rr$TGrouping)))
stat <- rr %>% filter(budget %in% c(10,30,50)) %>% 
  dplyr::group_by(run,feature,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T) ) %>% ungroup() %>% 
  dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
  # Average per run
  dplyr::group_by(feature,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T)) %>% ungroup() %>% 
  dplyr::mutate(gr = ifelse(feature %in% c('agb_carbon','water'),feature,'Biodiversity') )

# Shortfall for carbon and water globally
stat %>% filter(gr %in% c('agb_carbon','water')) %>% mutate(shortfall = 1 - shortfall)
stat %>% filter(gr %notin% c('AGB_carbon','Water')) %>% 
  dplyr::group_by(gr,budget) %>% 
  dplyr::summarise(avg_shortfall = mean(shortfall,na.rm = TRUE),
                   target_reached = sum(shortfall==0,na.rm = TRUE) / n() ) 

# What we get for 10% and 30%
rr_full %>% filter(budget %in% c(10,30,50))

# Amount of land for biodiversity only
rr_full %>% filter(type == 'Biodiversity only')

rr_full %>% filter(type == 'Biodiversity, Carbon and Water')

# Shortfall and targets reached of currently protected sites
rr %>% dplyr::filter(budget == '15') %>% 
  dplyr::group_by(run,feature,kingdom,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T)) %>% ungroup() %>% 
  group_by(budget) %>% 
  summarise(mean_shortfall = 1- mean(shortfall,na.rm = T),
            target_reached = sum(shortfall == 0,na.rm = T) / n() )

# Plants that have target reached at given budget
x = #rr %>% dplyr::filter(budget == '15',kingdom == "PLANTAE") %>%
  rr3 %>% dplyr::filter(budget == '10',kingdom == "PLANTAE") %>% # From bio+carbon+water
  dplyr::group_by(run,feature,kingdom,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T)) %>% ungroup() %>% 
  dplyr::group_by(feature,kingdom,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T)) %>% ungroup()
x %>% summarise(prop_plants = sum(shortfall==0,na.rm = TRUE) / n())

# Shortfall of carbon and water at 30 and 50%
rr %>% dplyr::filter(budget %in% c('10','30','50'), data %in% c('Carbon','Water')) %>% dplyr::group_by(budget,data) %>% dplyr::summarise(shortfall = 1-mean(shortfall))
# And improvement in conservation status
rr_full %>% dplyr::filter(budget %in% c('30','50')) %>% filter(type == 'Biodiversity, Carbon and Water')

# Summary per overall assets
rr %>% dplyr::mutate(gr = ifelse(feature %in% c('AGB_carbon','Water'),feature,'Biodiversity') ) %>% 
  dplyr::filter(budget %in% c('30','50')) %>% 
  dplyr::group_by(run,gr,feature,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T) ) %>% ungroup() %>% 
  # Average per feature
  dplyr::group_by(gr,feature,budget) %>% 
  dplyr::summarise(shortfall = mean(shortfall,na.rm = T) ) %>% ungroup() %>% 
  # Average per group
  dplyr::group_by(gr, budget) %>% dplyr::summarise(shortfall = mean(shortfall,na.rm=T))

# --- #
if(!repr_id){
  # Full comparison with all species
  rr_full <- bind_rows(
      rr %>% dplyr::group_by(run,feature,budget) %>% 
      dplyr::summarise(target_reached = sum(target_reached, na.rm = TRUE) / n() ) %>% ungroup() %>% 
      dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
      # Average + sd per run
      dplyr::group_by(budget) %>% 
      dplyr::summarise(target_reached_avg = mean(target_reached*100, na.rm = TRUE),
                       target_reached_ymin = mean(target_reached*100, na.rm = TRUE) - sd(target_reached*100, na.rm = TRUE),
                       target_reached_ymax = mean(target_reached*100, na.rm = TRUE) + sd(target_reached*100, na.rm = TRUE)) %>% 
      dplyr::mutate(type = 'Biodiversity only (representative)'),
    readRDS("temporary_SIallspecies_reprsave.rds") %>% 
      dplyr::group_by(run,feature,budget) %>% 
      dplyr::summarise(target_reached = sum(target_reached, na.rm = TRUE) / n() ) %>% ungroup() %>% 
      dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
      # Average + sd per run
      dplyr::group_by(budget) %>% 
      dplyr::summarise(target_reached_avg = mean(target_reached*100, na.rm = TRUE),
                       target_reached_ymin = mean(target_reached*100, na.rm = TRUE) - sd(target_reached*100, na.rm = TRUE),
                       target_reached_ymax = mean(target_reached*100, na.rm = TRUE) + sd(target_reached*100, na.rm = TRUE)) %>% 
      dplyr::mutate(type = 'Biodiversity only (all data)')
  )
}
# -------------------------- #
# Plotting below 
# -------------------------- #

# The plot for overall
g1 <- rr_full %>% 
  # ----- #
  #dplyr::filter(budget <= 60) %>% 
  ggplot(., aes(x = budget, y = target_reached_avg, group = type,colour = type)) + 
  theme_few(base_size = 20) +
#  geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
#  geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax), fill = 'grey80', alpha = .7) +
  geom_line(size = 2) +
  # Protected area indication
  geom_point(data = rr_full %>% filter(budget == 15), colour = 'red', alpha = .2, size = 3,stroke = 2) +
  scale_colour_manual(values = c("#00A087","#B09C85","#72bcd4","grey30"),guide = guide_legend("",override.aes = aes(fill = NA))) + #scale_colour_manual(values = 'black') + 
#  scale_linetype_manual(values = c('solid','dotted','longdash','twodash'),guide = guide_legend("",override.aes = aes(fill = NA))) + 
    theme(legend.position = c(.51,.15),legend.key.width = unit(1,'in'),legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent') ) +
  guides(fill = 'none') +
  scale_x_continuous(breaks = pretty_breaks(5),limits = c(as.numeric(levels(rr$budget)[1]),100)) +
  scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
  ylab("Targets reached (% of all species)") + xlab("Area managed for conservation\n(% of land area)") +
  theme(plot.background = element_blank())
g1
if(pa_fname == ""){
  ggsave(paste0(figure_path,"/Figure1_budgetrepresentation_overall.png"),plot =g1,dpi = 400)
} else {
  ggsave( paste0(figure_path,"/SIFigure1_budgetrepresentation_overall_withPA.png"),plot = g1,dpi = 400)
}

cols.tgrouping <- c(Plants = "#144720",Reptiles = "#871e7e",
                    Amphibians = "#1c98ba", Birds = "#b89216", Mammals = "#940303",
                    Carbon = 'grey10', Water = 'blue',
                    Overall = "#000000")

# Summarise per feature and group
rr_full_tgrouping <- rr %>% dplyr::group_by(TGrouping,feature,run,budget) %>% 
  dplyr::summarise(target_reached = sum(target_reached, na.rm = TRUE) / n() ) %>% ungroup() %>% 
  dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
  # Average + sd per run
  dplyr::group_by(TGrouping,feature,budget) %>% 
  dplyr::summarise(target_reached = mean(target_reached, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::group_by(TGrouping, budget) %>% 
  dplyr::summarise(target_reached_avg = mean(target_reached*100, na.rm = TRUE)
  ) 

# The plot for taxonomic groups
g2 <- rr_full_tgrouping %>% 
  # ----- #
  #dplyr::filter(budget <= 60) %>% 
  ggplot(., aes(x = budget, y = target_reached_avg, color = TGrouping, fill = TGrouping)) +
  theme_few(base_size = 20) +
#  geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
#  geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax), colour = NA, alpha = .4) + 
  geom_line(size = 2,linejoin = "round",lineend = "butt") +
  # Protected area indication
  geom_point(data = rr_full_tgrouping %>% filter(budget == 15), colour = 'red', alpha = .2, size = 3, stroke = 2) +
  scale_colour_manual(values = cols.tgrouping,guide = guide_legend("")) + theme(legend.position = c(.72,.15),legend.background = element_rect(fill = 'transparent') ) +
  scale_fill_manual(values = cols.tgrouping,guide = guide_legend("")) + guides(fill = 'none') +
  #scale_color_brewer(palette = "Set1",guide = guide_legend("")) + theme(legend.position = c(.85,.25)) +
  # Ablines
  #geom_vline(xintercept = 14,linetype = 'dotted',size=1) + annotate("label",x =17, y = 100, label = "2019") +
  #geom_vline(xintercept = 50,linetype = 'dotted',size=1) + annotate("label",x =53, y = 90, label = "Half-Earth") +
  scale_x_continuous(breaks = pretty_breaks(5),limits = c(as.numeric(levels(rr$budget)[1]),100)) +
  scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
  ylab("Targets reached (% of all assets)") + xlab("Areas managed for conservation\n(% of land area)") +
  theme(plot.background = element_blank())
g2
if(pa_fname == ""){
  ggsave(paste0(figure_path,"/Figure1_budgetrepresentation_tgrouping.png"),plot =g2,dpi = 400)
} else {
  ggsave(paste0(figure_path,"/SIFigure1_budgetrepresentation_tgrouping_withPA.png"),plot =g2,dpi = 400)
}

cols.category.threatened <- c(threatened = "#85052e",`non-threatened` = "#524f49")

# Overall
rr_full_threatened <- rr %>% 
  dplyr::filter(kingdom == "ANIMALIA") %>%
  dplyr::group_by(category.threatened,run,feature,budget) %>% 
  dplyr::summarise(target_reached = sum(target_reached, na.rm = TRUE) / n() ) %>% ungroup() %>% 
  dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
  # Average + sd per run
  dplyr::group_by(category.threatened,feature,budget) %>% 
  dplyr::summarise(target_reached = mean(target_reached, na.rm = TRUE)) %>% ungroup() %>% 
  dplyr::group_by(category.threatened, budget) %>% 
  dplyr::summarise(target_reached_avg = mean(target_reached*100, na.rm = TRUE)
#                   target_reached_ymin = mean(target_reached*100) - sd(target_reached*100),
#                   target_reached_ymax = mean(target_reached*100) + sd(target_reached*100)
                   ) 

# The plot for Threatened species
g3 <- rr_full_threatened %>% 
  # ----- #
  #dplyr::filter(budget <= 60) %>% 
  ggplot(., aes(x = budget, y = target_reached_avg, color = category.threatened,
                fill = category.threatened)) + 
  theme_few(base_size = 20) +
 # geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
 # geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax, fill = category.threatened), colour = NA, alpha = .4) + 
  geom_line(size = 2,linejoin = "round",lineend = "butt") +
  # Protected area indication
  geom_point(data = rr_full_threatened %>% filter(budget == 15), colour = 'red', alpha = .2, size = 3, stroke = 2) +
  scale_colour_manual(values = cols.category.threatened,guide = guide_legend("")) + theme(legend.position = c(.72,.15),legend.background = element_rect(fill = 'transparent') ) +
  scale_fill_manual(values = cols.category.threatened,guide = guide_legend("")) + guides(fill = 'none') +
  scale_x_continuous(breaks = pretty_breaks(5),limits = c(as.numeric(levels(rr$budget)[1]),100)) +
  scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
  ylab("Targets reached (% of all species)") + xlab("Areas managed for conservation\n(% of land area)") +
  theme(plot.background = element_blank())
g3
if(pa_fname == ""){
  ggsave(paste0(figure_path,"/Figure1_budgetrepresentation_threatened.png"),plot = g3,dpi = 400)
} else {
  ggsave(paste0(figure_path,"/SIFigure1_budgetrepresentation_threatened_withPA.png"),plot = g3,dpi = 400)
}
# Composite and export # 
library(cowplot)
remYaxis <- theme(axis.line.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())
pg <- plot_grid(g1,g2 + remYaxis,g3 + remYaxis ,nrow = 1,labels = 'auto',align = 'hv',axis = 'l', label_size = 24,label_fontface = 'bold')

if(pa_fname == ""){
  cowplot::ggsave2(filename = paste0(figure_path,"/Figure1_budgetrepresentation.png"),width = 18, height = 9,
                 plot = pg,dpi = 400)
} else {
  # Separate for PA lockedin
  cowplot::ggsave2(filename = paste0(figure_path,"/SIFigure1_budgetrepresentation_withPA.png"),width = 18, height = 9,
                   plot = pg,dpi = 400)
}

# ------------------------- #
# Assemble joint plot with protected areas locked in
out <- bind_rows(rr_full_forexport,rr_full_tgrouping,rr_full_threatened)
if(pa_fname == ""){ out1 <- out %>% mutate(Lockin_currentlyprotected = FALSE);rm(out) }
if(pa_fname != ""){ out2 <- out %>% mutate(Lockin_currentlyprotected = TRUE);rm(out) }
if(pa_fname == ""){ g1_woPA <- g1;g2_woPA <- g2; g3_woPA <- g3 }

# ------------------- #
write_csv(bind_rows(out1,out2),paste0("figures/Figure3_Extendenddata",split_id,".csv"))
stopifnot(pa_fname != "")

remYaxis <- theme(axis.line.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())

if(split_id != ""){
fig3 <- plot_grid( 
    # Overall
    g1_woPA + labs(x = "",title = "Optimal",tag = "a") + guides(colour = 'none',linetype = 'none') + theme(plot.title = element_text(hjust = .5)),
    g1 + labs(x = "", title = "Current protection",tag = "b") + remYaxis + guides(linetype = 'none',colour = 'none') + theme(plot.title = element_text(hjust = .5)),
    # Taxonomic groups
    g2_woPA + labs(x="",tag = "c") + guides(colour = 'none',linetype = 'none')  + theme(plot.title = element_text(hjust = .5)),
    g2 + labs(x="",tag = "d") + remYaxis + guides(linetype = 'none',colour = 'none'),
    # Threatened species
    g3_woPA + labs(tag = "e") + guides(colour = 'none',linetype = 'none')  + theme(plot.title = element_text(hjust = .5)),
    g3 + labs(tag = "f") +  remYaxis + guides(linetype = 'none',colour = 'none'),
    nrow = 3,align = 'hv',axis = 'l',label_x = .05, label_size = 20,label_fontface = 'bold')
  
cowplot::ggsave2(filename = paste0(figure_path,"/Figure3.png"),width = 14, height = 15,
                 plot = fig3,dpi = 400)
cowplot::ggsave2(filename = paste0(figure_path,"/Figure3.svg"),width = 14, height = 15,
                 plot = fig3,dpi = 400)
} else {
  fig3 <- plot_grid( 
    # Overall
    g1_woPA + labs(x = "",title = "Optimal",tag = "b") + guides(colour = 'none',linetype = 'none') + theme(plot.title = element_text(hjust = .5)),
    g1 + labs(x = "", title = "Current protection",tag = "c") + remYaxis + guides(linetype = 'none',colour = 'none') + theme(plot.title = element_text(hjust = .5)),
    # Taxonomic groups
    g2_woPA + labs(x="",tag = "d") + guides(colour = 'none',linetype = 'none')  + theme(plot.title = element_text(hjust = .5)),
    g2 + labs(x="",tag = "e") + remYaxis + guides(linetype = 'none',colour = 'none'),
    # Threatened species
    g3_woPA + labs(tag = "f") + guides(colour = 'none',linetype = 'none')  + theme(plot.title = element_text(hjust = .5)),
    g3 + labs(tag = "g") +  remYaxis + guides(linetype = 'none',colour = 'none'),
    nrow = 3,align = 'hv',axis = 'l',label_x = .05, label_size = 20,label_fontface = 'bold')
  
  cowplot::ggsave2(filename = paste0(figure_path,"/SIFigure10_nobiomesplit.png"),width = 14, height = 15,
                   plot = fig3,dpi = 400)
  cowplot::ggsave2(filename = paste0(figure_path,"/SIFigure10_nobiomesplit.svg"),width = 14, height = 15,
                   plot = fig3,dpi = 400)
}

# Save the legends separately and add them manually to the   plot
pl <- plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(g1 + theme(legend.key.width = unit(1,'in'))
                               )),
  ggpubr::as_ggplot(ggpubr::get_legend(g2 + theme(legend.key.width = unit(1,'in'))# +
#                                 guides(colour = guide_legend(title = "",label.position = 'top',direction = 'horizontal',
#                                                                label.hjust = 0.5,ncol=1,byrow = TRUE))
  )),
ggpubr::as_ggplot(ggpubr::get_legend(g3 + theme(legend.key.width = unit(1,'in'))# +
  #                               guides(colour = guide_legend(title = "",label.position = 'top',direction = 'horizontal',
  #                                                            label.hjust = 0.5,ncol=1,byrow = TRUE))
    ))
)
cowplot::ggsave2(filename = paste0(figure_path,"/Figure3_legends.svg"),width = 18, height = 15,
                 plot = pl,dpi = 400)

# ---------------------------- #
#### [Figure] Differing weights for carbon and water ####
# Idea:
# Run varying weights for carbon and water at 50km for a representative set
# Make a ternary plot with biodiversity, carbon and water and each axis
library(fst)
library(tidyverse)
library(ggplot2)
library(scales)
library(scico)
library(pbapply)
library(assertthat)

# Results path to the varying carbon_water weights
results_path = 'results/50km_esh_carbonwaterweights' # 
# List of files
ll_results <- list.files(results_path,full.names = T)
ll_results <- ll_results[has_extension(ll_results,'fst')]
ll_results <- ll_results[str_detect(ll_results,'securitysave',negate = T)] # Remove the securitysaves

# --- #
# Preprocessing function - Calculates shortfall and extracts the weight from filename
format_fname <- function(f){
  #print(basename(f))
  target_carbonweight = as.numeric( str_remove(str_split(tools::file_path_sans_ext( basename(f) ),"__",simplify = T)[,4],'carbweight') )
  target_waterweight = as.numeric( str_remove(str_split(tools::file_path_sans_ext( basename(f) ),"__",simplify = T)[,5],'waterweight') )
  target_budget = as.numeric( str_remove(str_split(tools::file_path_sans_ext( basename(f) ),"_",simplify = T)[,14],'perc') )
  # Load in File
  df <- read_fst(f) %>% 
    dplyr::filter(feature_abundance_pu > 0) %>% # No features with zero abundance
    dplyr::mutate(biocarbongroup = c(rep('Biodiversity',n()-2),'Carbon','Water')) %>% # Carbon and Water come last
    # Calculate the shortfall for each target / feature
    mutate(
      short_pu = pmax(absolute_target - absolute_held, 0)
    ) %>% 
    # express the shortfalls as a percentage of the targets
    mutate(
      shortfall = short_pu / absolute_target
    )
  if(str_detect(f,'_biome.id_')){ 
    df$feature <- str_split(df$feature,"__",simplify = T)[,1] 
#    df <- df %>% group_by(biocarbongroup,feature) %>% 
#      dplyr::summarise(absolute_held = mean(absolute_held),
#                       absolute_target = mean(absolute_target),
#                       shortfall = mean(shortfall)
#      )
    out <- df %>% group_by(biocarbongroup, feature) %>% 
      dplyr::summarise(shortfall = mean(shortfall,na.rm = TRUE) ) %>% ungroup() %>% 
      dplyr::group_by(biocarbongroup) %>% 
      dplyr::summarise(avg_shortfall = mean(shortfall,na.rm = TRUE),
                       target_reached = sum(shortfall==0,na.rm = TRUE) / n() ) %>% 
      dplyr::mutate(carbonweight = target_carbonweight,
                    waterweight = target_waterweight, budget = target_budget)      
  } else {
    # Now average and calculate
    out <- df %>% dplyr::group_by(biocarbongroup) %>% 
      dplyr::summarise(target_reached = length(which( absolute_held >= absolute_target )) / n() ,
                       avg_shortfall = mean(shortfall,na.rm = TRUE)) %>% ungroup() %>% 
      dplyr::mutate(carbonweight = target_carbonweight,
                    waterweight = target_waterweight, budget = target_budget)
  }
  return(out)
}

# Now do this for each file
xx <- pblapply(ll_results, function(x) {format_fname(x)} )
df <- do.call(rbind,xx)
df <- full_join(
  df %>% dplyr::filter(biocarbongroup == 'Biodiversity') %>% dplyr::select(-target_reached) %>%
    tidyr::spread(biocarbongroup,avg_shortfall),
  df %>% dplyr::select(-target_reached) %>%
    dplyr::filter(biocarbongroup %in% c("Carbon","Water")) %>%
    tidyr::spread(biocarbongroup,avg_shortfall),
  by = c("budget","carbonweight","waterweight")
)

unique(df$carbonweight)

# Reference / equal weight based on number of included assets
weight_equal <- 89618

# Remove those from xx that are smaller than equal_weight
df <- df %>% filter(waterweight <= weight_equal) %>% filter(carbonweight <= weight_equal)

# ---- #
# 2D plots for 30 and 50% budget globally and varying weights
# Show surface area? 
library(ggrepel)
xx <- df %>% 
  dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% #%>% dplyr::filter(budget %in% seq(10,70,10))
  #dplyr::mutate(Biodiversity = abs(Biodiversity - 100)) %>% 
  # Add overall score
  rowwise() %>%  dplyr::mutate(score = sum(Biodiversity,Carbon,Water))
# Stats for ranges
range(xx$carbonweight)

# Biodiv only
xx %>% dplyr::filter(budget == 10, carbonweight == weight_equal,waterweight == 1)
xx %>% dplyr::filter(budget == 10, carbonweight == 1,waterweight == weight_equal)

xx %>% dplyr::filter(budget == 30,carbonweight == weight_equal,waterweight == 1)
xx %>% dplyr::filter(budget == 30,carbonweight == 1,waterweight == weight_equal)

xx %>% #slice(11:n()) %>% 
  dplyr::filter(budget %in% c(10,30)) %>% group_by(budget) %>% 
  summarise(
    min_carbon = min(Carbon),  max_carbon = max(Carbon),
    min_water = min(Water),  max_water = max(Water)
)
xx %>% dplyr::filter(budget == 30, carbonweight == weight_equal,waterweight == weight_equal)

# Largest aggregated score
xx %>% dplyr::filter(carbonweight == weight_equal,waterweight == weight_equal) %>% arrange(desc(score))

# --- #
xx_bioonly <- xx %>% dplyr::slice(1:10)
xx <- xx %>% dplyr::slice(11:n())
write_csv(xx,paste0('figures/Figure2_Extendenddata',split_id,'.csv'))

# Weighting as size
xx = xx %>% rowwise() %>%  dplyr::mutate(weight_diff = mean(c(carbonweight,waterweight)) ) %>% ungroup()
# Convert into factor/ranks
xx$weight_diff <- dense_rank(xx$weight_diff)
# Custom break labels
frac_lab <- data.frame( b = sort(unique(xx$carbonweight)),
                        r = seq(1,10),
                        l = c("None",
#                              "~1/10",
                              "~2/10",
                              "~3/10",
                              "~4/10",
                              "~5/10",
                              "~6/10",
                              "~7/10",
                              "~8/10",
                              "~9/10",
                              "Equal"
                        )
)

# Highlighted points
xx$highlight <- NA
xx$highlight[intersect( which(xx$carbonweight==weight_equal) , which( xx$waterweight==1))] <- 1
xx$highlight[intersect( which(xx$carbonweight==1) , which( xx$waterweight==weight_equal))] <- 1
xx$highlight[which(xx$budget %notin% c(10,30))] <- NA
xx %>% dplyr::filter(highlight == 1)

gcop <- ggplot(xx, aes(x = Carbon, y = Water, colour = Biodiversity, size = weight_diff, group = weight_diff)) +
    theme_bw(base_size = 20,base_family = 'Arial') +
    geom_point(data = subset(xx, highlight == 1),aes(size = weight_diff+0.1), colour = 'black', stroke = 1.25) +
    # Biodiversity only
    geom_point(alpha = .85,fill = NA) +
    # # Labels for specific points
    geom_text_repel(
      data          = subset(xx, highlight == 1),
      label         = c('(c)','(c)','(b)','(b)'),
      colour       = "black", size = 5,
      nudge_x      = 7, xlim = c(NA,80),
      nudge_y       = 5,
      segment.size  = 0.35,
      segment.color = "black",
      direction     = "both",
      point.padding = NA,
      show.legend = FALSE
    ) +
    scale_size(#breaks = c(1,3,5,7,11),labels = frac_lab$l[c(1,3,5,7,11)],
               breaks = c(1,2,5,7,10),labels = frac_lab$l[c(1,2,5,7,10)],
               #breaks = frac_lab$r, labels = frac_lab$l,
               range = c(1,7), guide = guide_legend(title = "Weighting",title.theme = element_text(size = 14),label.position = 'bottom',
                                                                             direction = 'horizontal',title.position = 'top',
                                                                             label.theme = element_text(size = 12)) ) +
    # Colours
    scale_colour_scico(10,palette = 'roma',breaks = c(65,75,85,95,100),guide = guide_colourbar(title = 'Improvement in conservation status (%)',#"Biodiversity shortfall reduced (%)",
                                                                                               title.position = 'top',
                                                                   direction = 'horizontal',barwidth = grid::unit(3,'in'),
                                                                   title.theme = element_text(size = 14),label.theme = element_text(size = 12), ticks = FALSE)) +
    theme(legend.position = c(.75,.20), legend.background = element_rect(fill = 'transparent'),legend.key =  element_rect(fill = 'transparent',colour = NA)) +
    # Arrow
    labs(x = "Total carbon conserved (%)",
         y = "Total water conserved (%)"
         )
gcop
ggsave(filename = paste0(figure_path,"/","Figure2_varyingweights.png"),plot = gcop,width = 9, height = 9,dpi = 400)
ggsave(filename = paste0(figure_path,"/","Figure2_varyingweights.svg"),plot = gcop,width = 9, height = 9,dpi = 400)

# Split the results into varying carbon and water weights
df.carbon <- bind_rows( df %>% dplyr::filter(carbonweight < weight_equal), df %>% dplyr::filter(carbonweight == weight_equal, waterweight == weight_equal) )
df.water <- bind_rows( df %>% dplyr::filter(waterweight < weight_equal), df %>% dplyr::filter(carbonweight == weight_equal, waterweight == weight_equal) )

df.carbon <- df.carbon %>% dplyr::mutate_at(c('Carbon','Water'),  function(x) { abs(x- 1)*100 } ) %>% 
  dplyr::mutate(Biodiversity = Biodiversity * 100)
df.water <- df.water %>% dplyr::mutate_at(c('Carbon','Water'),  function(x) { abs(x- 1)*100 } ) %>% 
  dplyr::mutate(Biodiversity = Biodiversity * 100)


# ---- #
# Double axes plot for each budget up to 70%

# Carbon weight only
x1 <- xx %>% slice(1:90,180:190) %>% filter(budget<= 70)
# Water Weights only
x2 <- xx %>% slice(91:190) %>% filter(budget<= 70)

biPlot <- function(x, y, z, title, coly = 'black', colz = 'blue'){
  # add extra space to right margin of plot within frame
  par(mar=c(5, 4, 4, 6) + 0.01)
  
  assert_that(length(x) == length(y), length(y) == length(z))
  # Axis limit
  ml = round(ifelse(min(y,z)==0,0, min(y,z)-10 ),-1)
  ul = round(ifelse(max(y,z)>=100,100, max(y,z)+10 ),-1)
  ulb = round(ifelse(max(x)>=100,100, max(x)+10 ),-1)
  
  # Plot first set of data and draw its axis
  plot(x = x, y = y, pch=16, axes=FALSE, ylim=c(ml,ul), xlab="", ylab="", 
       type="b",col="black", main= title)
  axis(2, ylim = c(ml,ul), col="black",las=1)  ## las=1 makes horizontal labels
  mtext("Total carbon conserved (%)",side=2,line=2.5)
  box()
  
  # Allow a second plot on the same graph
  par(new=TRUE)
  
  # Plot the second plot and put axis scale on right
  plot(x = x, y = z, pch=15,  xlab="", ylab="", ylim=c(ml,ul), 
       axes=FALSE, type="b", col="blue")
  # a little farther out (line=4) to make room for labels
  mtext("Total water conserved (%)",side=4,col="blue",line=3) 
  axis(4, ylim = c(ml,ul), outer = FALSE, col="blue",col.axis="blue",las=1)
  
  # Draw the biodiversity axes
  axis(1, xlim = pretty(c(60, 100),10) )
  mtext("Biodiversity shortfall reduced (%)",side=1,col="black",line=2.5)
}

# Now build plots
png('figures/SIFigure_ShortfallWeights.png',width = 3000,height = 3600,units = "px",res = 400)
par(mfrow = c(3,2))
for(b in c(10,30,50)){
  
  sub1 <- subset(x1,budget == b)
  sub2 <- subset(x2,budget == b)
  
  biPlot(x = sub1$Biodiversity, y = sub1$Carbon, z = sub1$Water,
         title = paste0("Varying carbon weight\n for ",b, "% of land area" ))
  biPlot(x = sub2$Biodiversity, y = sub2$Carbon, z = sub2$Water,
         title = paste0("Varying water weight\n for ",b, "% of land area"))
}
dev.off()
      

#### [Figure] Maps and spiderplots for differential weights of all assets ####
# Idea: 
# Show 30% and 50% for all variantes
# Biodiversity full, Carbon and water 1
# Biodiversity full, Carbon full, water 1
# Biodiversity full, Carbon 1, water full
library(raster)
library(ggradar)
library(cowplot)
library(ggpmisc)

# The equal weighting
weight_equal <- 89618 #198998

stopifnot( assert_that( exists('format_fname'),exists('ll_results') ) )  # Load from above

xx <- pblapply(ll_results[union(grep('10perc',ll_results),grep('30perc',ll_results))], function(x) {format_fname(x)} )
df <- do.call(rbind,xx);rm(xx)
df <- full_join(
  df %>% dplyr::filter(biocarbongroup == 'Biodiversity') %>% dplyr::select(-avg_shortfall) %>% 
    tidyr::spread(biocarbongroup,target_reached),
  df %>% dplyr::select(-target_reached) %>% 
    dplyr::filter(biocarbongroup %in% c("Carbon","Water")) %>% 
    tidyr::spread(biocarbongroup,avg_shortfall),
  by = c("budget","carbonweight","waterweight")
)

# Load representation files and make 3d plot
#biodiv_only <- df %>% dplyr::filter(budget %in% c(10,30)) %>% 
#  dplyr::filter(carbonweight == 1, waterweight == 1)
carbon_only <- df %>% dplyr::filter(budget %in% c(10,30)) %>% 
  dplyr::filter(carbonweight == weight_equal, waterweight == 1)
water_only <- df %>% dplyr::filter(budget %in% c(10,30)) %>% 
  dplyr::filter(waterweight == weight_equal, carbonweight == 1)
stopifnot(nrow(carbon_only)>1,nrow(water_only)>1 )

# Colours
cols = c( scico(10,palette = 'roma')[1],scico(10,palette = 'roma')[3] )

# # Make plot for biodiversity only
# gr1 <- biodiv_only %>% dplyr::select(budget,Biodiversity:Water) %>% 
#   dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% 
#   dplyr::mutate(Biodiversity = abs(Biodiversity - 100)) %>% 
#   ggradar(.,grid.min = 0,grid.mid = 50, grid.max = 100,
#           axis.label.size = 6, grid.label.size = 4, legend.text.size = 10,
#           base.size = 4
#           ) +
#   theme_void(base_size = 8) +
#   #  labs(
#   #    x = "Species at risk of extinction (%)",
#   #    y = "Total carbon conserved (%)",
#   #    z = "Total water conserved (%)" 
#   #  )
#   scale_colour_manual(values = cols) +
#   guides(colour = 'none')
# gr1

gr2 <- carbon_only %>% dplyr::select(budget,Biodiversity:Water) %>% 
  dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% 
  dplyr::mutate(Biodiversity = abs(Biodiversity - 100)) %>% 
  ggradar(.,
          grid.min = 0,grid.mid = 50, grid.max = 100,
          axis.label.size = 6, grid.label.size = 7, gridline.label.offset = 0,legend.text.size = 10,
          base.size = 20, group.point.size = 4,group.line.width = 1.05
          #          axis.labels = c('Biodiversity','Carbon','Water'), plot.title = ''
  ) +
  theme_void(base_size = 8) +
  theme(legend.title = element_text(size = 12, angle = 90),legend.text = element_text(size = 12),
        legend.position = 'left') + guides(colour = guide_legend(title = "Total land area (%)",nrow = 1)) +
  #  labs(
  #    x = "Species at risk of extinction (%)",
  #    y = "Total carbon conserved (%)",
  #    z = "Total water conserved (%)" 
  #  )
  scale_colour_manual(values = cols) +
  guides(colour = 'none')
gr2
ggsave(filename = paste0(figure_path,"/","Figure2_legend_carbononly.png"),plot = gr2, width = 5, height = 5, dpi = 400,bg = 'transparent')
ggsave(filename = paste0(figure_path,"/","Figure2_legend_carbononly.svg"),plot = gr2, width = 5, height = 5, dpi = 400,bg = 'transparent')

gr3 <- water_only %>% dplyr::select(budget,Biodiversity:Water) %>% 
  dplyr::mutate_at(c('Biodiversity','Carbon','Water'),  function(x) { abs( (x * 100) - 100 ) } ) %>% 
  dplyr::mutate(Biodiversity = abs(Biodiversity - 100)) %>% 
  ggradar(.,
          grid.min = 0,grid.mid = 50, grid.max = 100,
          axis.label.size = 6, grid.label.size = 7, gridline.label.offset = 0,legend.text.size = 10,
          base.size = 20, group.point.size = 4,group.line.width = 1.05
          #          axis.labels = c('Biodiversity','Carbon','Water'), plot.title = ''
  ) +
  theme_void(base_size = 8) +
  #  labs(
  #    x = "Species at risk of extinction (%)",
  #    y = "Total carbon conserved (%)",
  #    z = "Total water conserved (%)" 
  #  )
  scale_colour_manual(values = cols) +
  guides(colour = 'none')
gr3
ggsave(filename = paste0(figure_path,"/","Figure2_legend_wateronly.png"),plot = gr3, width = 5, height = 5, dpi = 400,bg = 'transparent')
ggsave(filename = paste0(figure_path,"/","Figure2_legend_wateronly.svg"),plot = gr3, width = 5, height = 5, dpi = 400,bg = 'transparent')

# List of tif files
results_path = 'results/50km_esh_carbonwaterweights'
ll_results <- list.files(results_path,full.names = T) # List of files
ll_results <- ll_results[has_extension(ll_results,'tif')]
ll_results <- ll_results[union(grep('10perc',ll_results),grep('30perc',ll_results))] # Only budgets 30 and 50%
# Grep the respective maps
ras_biodiv_only <- raster::stack( ll_results[intersect(grep("carbweight1_",ll_results), grep("waterweight1_",ll_results))] )
ras_carbon_only <- raster::stack( ll_results[intersect(grep(paste0("carbweight",weight_equal),ll_results), grep("waterweight1_",ll_results))] )
ras_water_only <- raster::stack( ll_results[intersect(grep("carbweight1_",ll_results), grep(paste0("waterweight",weight_equal),ll_results))] )

stopifnot(
  assert_that(nlayers(ras_carbon_only)==2,
              nlayers(ras_water_only)==2)
          )
# Get land area for continent outline
land <- raster('data/globalgrid_mollweide_50km.tif')
land[land==0] <- NA;land[land>0] <- 0

# Make a ranked map of 30% and 50%. Set to binary, since differences are minor and not as importanct for the figure
ras_biodiv_only <- sum( reclassify(ras_biodiv_only,c(0,1,1) ),na.rm = TRUE )
ras_biodiv_only[ras_biodiv_only ==0] <- NA
ras_biodiv_only <- mask(sum(ras_biodiv_only,land,na.rm = T),land)
ras_carbon_only <- sum( reclassify(ras_carbon_only,c(0,1,1) ),na.rm = TRUE )
ras_carbon_only[ras_carbon_only ==0] <- NA
ras_carbon_only <- mask(sum(ras_carbon_only,land,na.rm = T),land)
ras_water_only <- sum( reclassify(ras_water_only,c(0,1,1) ),na.rm = TRUE )
ras_water_only[ras_water_only ==0] <- NA
ras_water_only <- mask(sum(ras_water_only,land,na.rm = T),land)
# Stack them
ras_stack = stack(ras_biodiv_only,ras_carbon_only,ras_water_only)
names(ras_stack) <- c('Biodiversity only','Carbon preference','Water preference')

# ----- #
no_background <- theme(plot.background = element_rect( fill =  'transparent'),panel.background = element_rect(fill = "transparent"),legend.background = element_rect(fill = 'transparent') )

cols = c( scico(10,palette = 'roma')[1],scico(10,palette = 'roma')[4],'#F2F2F2' )

# Biodiv only
(out_fname <- paste0(figure_path, "/","Figure2_solution_biodivonly", ".png" ))
g1 <- ggR(ras_stack, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE, forceCat = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_manual(values = rev(cols), na.value = 'transparent') + guides(fill = 'none') +
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = out_fname, plot = g1,width = 8, height = 8, dpi = 400,bg = 'transparent')
# Make map with inset plot
g1_inset <- ggdraw() +
  draw_plot(g1) +
  draw_plot(gr1 + no_background + panel_border(remove = TRUE),
            x = 0.1, y = .25, width = .25, height = .25)
cowplot::ggsave2(plot = g1_inset, filename = out_fname,width = 14, height = 8, dpi = 400)

# Carbon only
(out_fname <- paste0(figure_path, "/","Figure2_solution_carbononly", ".png" ))
g2 <- ggR(ras_stack, layer = "Carbon.preference", maxpixels = 1e10,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE, forceCat = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_manual(values = rev(cols), na.value = 'transparent') + guides(fill = 'none') +
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "")
ggsave(filename = out_fname, plot = g2,width = 8, height = 8, dpi = 400,bg = 'transparent')
# Make map with inset plot
g2_inset <- ggdraw() +
  draw_plot(g2) +
  draw_plot(gr2 + no_background + panel_border(remove = TRUE),
            x = 0.1, y = .25, width = .25, height = .25)
cowplot::ggsave2(plot = g2_inset, filename = out_fname,width = 14, height = 8, dpi = 400)

# Water only
(out_fname <- paste0(figure_path, "/","Figure2_solution_wateronly", ".png" ))
g3 <- ggR(ras_stack, layer = "Water.preference", maxpixels = 1e10,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE, forceCat = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_manual(values = rev(cols), na.value = 'transparent') + guides(fill = 'none') +
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = out_fname, plot = g3,width = 8, height = 8, dpi = 400,bg = 'transparent')
# Make map with inset plot
g3_inset <- ggdraw() +
  draw_plot(g3) +
  draw_plot(gr3 + no_background + panel_border(remove = TRUE),
            x = 0.1, y = .25, width = .25, height = .25)
cowplot::ggsave2(plot = g3_inset, filename = out_fname,width = 14, height = 8, dpi = 400)

#### [SI-Figure] Global ranked carbon and water assets only ####
carbon_agbc <- raster( paste0(carbon_path,"/","Carbon_AGBC_",target_resolution,"_",projection,".tif") )
carbon_bgbc <- raster( paste0(carbon_path,"/","Carbon_BGBC_",target_resolution,"_",projection,".tif") )
carbon_soc <- raster( paste0(carbon_path,"/","Carbon_SOC_",target_resolution,"_",projection,".tif") )
# Align them all
if(!compareRaster(carbon_agbc,carbon_soc,stopiffalse = F)){
  # First for SOC
  carbon_soc <- raster::crop(carbon_soc,  land)
  carbon_soc <- extend(carbon_soc, land)
  carbon_soc <- raster::resample(carbon_soc,land,method = "bilinear")
  # Now for AGBC
  carbon_agbc <- raster::crop(carbon_agbc,  land)
  carbon_agbc <- extend(carbon_agbc, land)
  carbon_agbc <- raster::resample(carbon_agbc,land,method = "bilinear")
}
# Security checks
stopifnot(
  assert_that(
    compareRaster(carbon_soc,carbon_agbc),
    ncell(land) == ncell(carbon_soc),
    compareRaster(land,carbon_soc)
  )
)

# Filename for 
(fname_carb <- paste0(results_path,"/Carbon_",target_resolution,"_",projection,".tif"))
# Aggregate both layers
carb <- stack(carbon_soc,carbon_agbc,carbon_bgbc)
carb <- sum(carb,na.rm = TRUE)
carb[carb == 0] <- NA # Overwrite with NA
carb <- raster::mask(carb, land)
writeGeoTiff(carb, fname = fname_carb,'FLT4S') # Save outputs and clear 

(fname_carb <- paste0(results_path,"/Carbon_normalized_ranked_",target_resolution,"_",projection,".tif"))
carb_n <- normalize_raster(carb)*100
carb_n <- abs( raster_to_ranks(carb_n,n = 100) - 100 )
carb_n <- raster::mask(carb_n, land.mask)
#plot(carb_n,col = rev( RColorBrewer::brewer.pal(10,"RdYlBu") ))
writeGeoTiff(carb_n, fname = fname_carb) # Save outputs and clear 

# (fname_carb <- paste0(results_path,"/Carbon_normalized_",target_resolution,"_",projection,".tif"))
# carb_n <- abs( (normalize_raster(carb)*100) - 100 )
# globalgrid_frac = raster( paste0("data/globalgrid_",projection,"_",target_resolution,".tif" )) / 1000
# globalgrid_frac[globalgrid_frac>0.5] <- 1
# globalgrid_frac[globalgrid_frac<=0.5] <- NA
# carb_n <- raster::mask(carb_n, globalgrid_frac)
# writeGeoTiff(carb_n, fname = fname_carb) # Save outputs and clear 

# Plot the carbon layer
g1 <- ggR(carb_n, layer = 1, maxpixels = 1e10,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_gradientn(colours = scico(10,palette = 'oleron'),na.value = 'transparent',#breaks = c(1,25,50,75,100),
                       guide = guide_colourbar(title = "Priority rank",title.theme = element_text(size = 14),
                                               title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                               direction = 'horizontal',ticks = F,reverse = T) ) + 
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = paste0('figures/SIFigure3_Carbon_only_normalized_',target_resolution,"_",projection,'.png'),
       plot = g1,width = 8, height = 8, dpi = 400)

# Also for water
water <- raster( paste0(water_path,"/","cleanwaternorm_pot_",target_resolution,"_",projection,".tif") )
# Security checks
stopifnot(
  assert_that(
    compareRaster(water,land)
  )
)
(fname_water <- paste0(results_path,"/Water_normalized_ranked_",target_resolution,"_",projection,".tif"))
water[water==0] <- NA # No water here
water_n <- normalize_raster(water)*100
water_n <- abs( raster_to_ranks(water_n,n = 100) - 100 )
water_n <- raster::mask(water_n, land.mask)
#plot(water_n,col = rev( RColorBrewer::brewer.pal(10,"RdYlBu") ))
writeGeoTiff(water_n, fname = fname_water) # Save outputs and clear 

g2 <- ggR(water_n, layer = 1, maxpixels = 1e10,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_gradientn(colours = scico(10,palette = 'devon'),na.value = 'transparent',breaks = c(1,25,50,75,100),
                       guide = guide_colourbar(title = "Priority rank",title.theme = element_text(size = 14),
                                               title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                               direction = 'horizontal',ticks = F,reverse = T) ) + 
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = paste0('figures/SIFigure3_Water_only_normalized_',target_resolution,"_",projection,'.png'),
       plot = g2,width = 8, height = 8, dpi = 400)


# Carbon prioritization with protected areas locked in
ll_carbon <- list.files('results/10km_esh_carbononly/',pattern = '*.tif',full.names = TRUE)

ras2 <- makeStack(ll_carbon)
#tifs <- grep(c('50perc'),ll.tif,invert = F,value = T)
ras2 <- reclassify(ras2, cbind(NA, NA, 0), right = FALSE)
mcp_ilp_hier <- mean(ras2, na.rm = T) # Average stack
# Mask with global land area mask
mcp_ilp_hier <- raster::mask(mcp_ilp_hier, land.mask)

mcp_ilp_hier_norm <- mcp_ilp_hier * 100 # precision of input layers
# -- #
# Remove 0 fractional estimates
mcp_ilp_hier_norm[mcp_ilp_hier_norm==0] <- NA
# Convert to ranks and invert
mcp_ilp_hier_norm <- abs( raster_to_ranks(mcp_ilp_hier_norm,n = 100,plot = FALSE) - 100 )
# Mask again with global land area mask 
mcp_ilp_hier_norm[is.na(mcp_ilp_hier_norm)] <- 100 # Add lest important areas as 100 so that proportions match again
mcp_ilp_hier_norm <-  raster::mask(mcp_ilp_hier_norm,land.mask)
NAvalue(mcp_ilp_hier_norm) <- -9999
plot(mcp_ilp_hier_norm,col = scico::scico(10,palette = 'roma'))

# Save outputs #
fname = 'results/10km_esh_ranked/minshort_carbon_withPA_10km_ranked.tif'
writeGeoTiff(mcp_ilp_hier_norm, fname = fname) # Save outputs and clear  


g <- ggR(mcp_ilp_hier_norm, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_gradientn(breaks = c(100,75,50,25,1),colours = scico(10,palette = 'bamako'),na.value = 'transparent',
                       guide = guide_colourbar(title = "Priority rank",title.theme = element_text(size = 14),
                                               title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                               direction = 'horizontal',ticks = F,reverse = T) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
ggsave(filename = 'figures/Figure_carbon_withPA_10km_ranked.png', plot = g,width = 8, height = 8, dpi = 400)

####  [SI-Figure] Resample finer resolution map to 50km and compare! ####
# Accucary comparison betwen low and high resolution layers
library(greenbrown);library(gdalUtils)
library(ggthemes)
r1 <- raster("results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked_majorityclipped.tif")
r2 <- raster("results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked_majorityclipped.tif")
# Calculate average at 50%
gdal_translate(src_dataset = r1@file@name,
               dst_dataset = '10km_biocarbonwaterranked_at50km.tif',
               ot = 'Float32',
               tr = res(r2),
               co = c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
               multi = TRUE,
               r = 'average'
)
# Load and align biome layer
biomes <- raster('/media/martin/data/constraints/biomes_10km_mollweide_modal.tif')
biomes <- alignRasters(biomes, r2, method = 'ngb', func = raster::modal)

r1 <- raster('10km_biocarbonwaterranked_at50km.tif')
r1 <- alignRasters(r1,r2,method = 'ngb',func = raster::modal)
assert_that(
  compareRaster(r1,r2), compareRaster(r1,biomes)
)
r1[is.na(r2)] <- NA

land <- raster('data/globalgrid_mollweide_50km.tif') / 1000

# Create a paired comparison of the maps next to each other
difference = (r2 - r1)
cellStats(difference,'mean');cellStats(difference,'sd')
# The plot
g1 <- ggR(difference, layer = 1, maxpixels = 1e10,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  cowplot::theme_nothing(font_size = 11,font_family = 'Arial') +
    scico::scale_fill_scico(palette = 'nuuk',na.value = 'transparent',
                             guide = guide_colourbar(title = "Absolute difference\n(50km-10km)",title.theme = element_text(size = 14),ticks = F,reverse = TRUE)) +
#  scale_fill_gradient2(na.value = 'transparent',high = 'grey50',low = muted("red"),mid = 'grey95',
#                       guide = guide_colourbar(title = "Absolute difference\nin asset rank",frame.colour = NULL,title.theme = element_text(size = 14),ticks = F,reverse = TRUE)) +
  theme(legend.position = c(.2,.3),legend.background = element_rect(fill = 'transparent',colour = NA)) +
  labs(x = "", y = "" )
g1

# Idea 2:
# Calculate the accucary at each hieracichal level and plot as barplot
df <- data.frame(biome = getValues(biomes),
                 land = getValues(land),
                 asset10km = getValues(r1),
                 asset50km = getValues(r2)
                 ) %>% drop_na()

# Get overall area as weight
# Mean absolute error
MAE <- function(obs,pred) { mean(  (obs - pred) )   }

# Now for each group, calculate the MSE
out <- df %>% dplyr::mutate(diff = asset50km - asset10km) %>% 
  dplyr::group_by(biome) %>% 
  dplyr::summarise(difference = sum(diff * land) / ( n() * sum(land) ) ) 
# The biomes code
biomes_code <- sf::read_sf('/media/martin/data/raw/Ecoregions2017/Ecoregions2017.shp') %>% 
  dplyr::select(BIOME_NUM,BIOME_NAME,COLOR_BIO) %>% sf::st_drop_geometry() %>% distinct() %>% 
  dplyr::slice(-15) %>% 
  dplyr::arrange(BIOME_NUM)
biomes_code$BIOME_NAME <- factor(biomes_code$BIOME_NAME,levels = biomes_code$BIOME_NAME) # Convert to factor using current order

# Left_join with calculation
out <- left_join(out, biomes_code, by = c('biome' = 'BIOME_NUM'))

# Make the plot
g2 <- ggplot(out, aes(x = BIOME_NAME, y = difference, fill = BIOME_NAME)) +
  theme_few(base_size = 16) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = out$COLOR_BIO,guide = guide_legend(title = "",ncol = 3)) +
  scale_y_continuous(breaks = pretty_breaks(5),expand = c(0,0) ) +
  scale_x_discrete(position = "top") + theme(axis.text.x.top = element_text(hjust = 0,angle = 45)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'bottom',legend.text = element_text(size = 10)) +
  labs(x = "", y = "Mean difference (50km to 10km)")
g2

library(cowplot);library(gridExtra)
pg1 <- grid.arrange(g1,g2)
ggsave(filename = paste0(figure_path,'/SIFigure_scalecomparison.png'), plot = pg1,
       width = 12, height = 10,dpi = 400 )

pg2 <- plot_grid(g1,g2,
                labels = 'auto',label_size = 18,rel_widths = c(1,.5),rel_heights = c(1,1),rows = 2)
cowplot::ggsave2(filename = paste0(figure_path,'/SIFigure_scalecomparison.png'),plot = pg2,
                   width = 12, height = 12,dpi = 400)


####  [SI-Figure] Alternative - Resolution comparison between 10 km and 50 km ####
# Compare side by side for 10 km and 50 km
library(patchwork)
library(cowplot)

bio_only_10km <- raster('results/10km_esh_ranked/minshort_speciestargets_esh10km_repruns10_ranked_majorityclipped.tif')
bio_only_50km <- raster('results/50km_esh_ranked/minshort_speciestargets_esh50km_repruns10_ranked_majorityclipped.tif')

biocarb_10km <- raster('results/10km_esh_ranked/minshort_speciestargets_carbon__esh10km_repruns10_ranked_majorityclipped.tif')
biocarb_50km <- raster('results/50km_esh_ranked/minshort_speciestargets_carbon__esh50km_repruns10_ranked_majorityclipped.tif')

biocarbwater_10km <- raster('results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked_majorityclipped.tif')
biocarbwater_50km <- raster('results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked_majorityclipped.tif')

# Make a function to plot them
makePlot <- function(ras){
  ggR(ras, layer = 1, maxpixels = 1e10,
           geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
    ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
    scale_fill_gradientn(colours = scico(10,palette = 'roma'),na.value = 'transparent',
                         guide = guide_colourbar(title = "Asset rank",title.theme = element_text(size = 14),ticks = F,reverse = TRUE) ) +
    theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
    labs(x = "", y = "" ) + guides(fill = 'none')
}
calcRepr <- function(rr){
  rr %>% dplyr::group_by(run,budget) %>% 
    dplyr::summarise(target_reached = sum(target_reached) / n() ) %>% ungroup() %>% 
    dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
    # Average + sd per run
    dplyr::group_by(budget) %>% 
    dplyr::summarise(target_reached_avg = mean(target_reached*100),
                     target_reached_ymin = mean(target_reached*100) - sd(target_reached*100),
                     target_reached_ymax = mean(target_reached*100) + sd(target_reached*100))
}
makeReprPlot <- function(rr){
  ggplot(rr, aes(x = budget, y = target_reached_avg, group = resolution, colour = resolution)) + 
    theme_few(base_size = 20) +
    #  geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
    geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax), fill = 'grey80', alpha = .7) +
    geom_line(size = 2) +
    #scale_linetype_manual(values = c('solid','longdash'),guide = guide_legend("",override.aes = aes(fill = NA))) + 
    scale_colour_manual(values = c("blue","orange"),guide = guide_legend("",override.aes = aes(fill = NA))) + 
    guides(fill = 'none') + 
    theme(legend.position = c(.51,.15),legend.key.width = unit(1.5,'in'),legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent',colour = NA) ) +
    scale_x_continuous(breaks = pretty_breaks(5),limits = c(10,100)) +
    scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
    ylab("Targets reached (% of all species)") + xlab("Area managed for conservation\n(% of land area)") +
    theme(plot.background = element_blank())
}

g1 <- makePlot(bio_only_10km)
g2 <- makePlot(bio_only_50km)

g4 <- makePlot(biocarb_10km)
g5 <- makePlot(biocarb_50km)

g7 <- makePlot(biocarbwater_10km)
g8 <- makePlot(biocarbwater_50km)

# ------------------------------------------- #
# Also insert the representation curves for each comparison
# Bio only
rr_b_10km <- read_rds('temporary_bioonly_10km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '10km')
rr_b_50km <- read_rds('temporary_bioonly_50km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '50km')
# Biocarbon
rr_bc_10km <- read_rds('temporary_biocarbon_10km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '10km')
rr_bc_50km <- read_rds('temporary_biocarbon_50km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '50km')
# Biocarbonwater
rr_bcw_10km <- read_rds('temporary_biocarbonwater_10km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '10km')
rr_bcw_50km <- read_rds('temporary_biocarbonwater_50km_reprsave.rds') %>% calcRepr() %>% mutate(resolution = '50km')

# The plot for accumulation curves
g3 <- makeReprPlot( bind_rows(rr_b_10km,rr_b_50km) )
g6 <- makeReprPlot( bind_rows(rr_bc_10km,rr_bc_50km) )
g9 <- makeReprPlot( bind_rows(rr_bcw_10km,rr_bcw_50km) )
# ----- #

# Assemble
pg <- plot_grid(g1, g2, g3, g4, g5, g6, 
                align = 'hv',ncol = 2,scale = 1.25)
cowplot::ggsave2(filename = paste0('figures/SIFigure9_comparison.png'), plot = pg,width = 14, height = 16, dpi = 400)
cowplot::ggsave2(filename = paste0('figures/SIFigure9_comparison_raw.svg'), plot = pg,width = 14, height = 16, dpi = 400)

pg <- plot_grid(g3 + xlab(''),g6 + xlab(''),g9,nrow = 3, align = 'hv')
cowplot::ggsave2(filename = paste0('figures/SIFigure9_comparison_reprcurves.png'), plot = pg,
                 width = 1000, height = 3300,units = 'px', dpi = 400)


# ------------------------------------------------------- #
#### [SI-Figure] Comparison of species weights for Threat,ED and EDGE ####
# Idea: Take the no weight raster as baseline and 10% (30%) from each species
# Overlay them all in a Venn diagram. Also consider   
library(RStoolbox);library(cowplot);library(gridExtra)

what = "10km"
pu_id_raster <- raster(str_glue('data/puid_{what}.tif'))
land <- raster(str_glue('data/globalgrid_mollweide_{what}.tif'))
ras_normal_comparison <- raster(str_glue('results/10km_esh/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh10km_10perc.tif'))
ras_threat_comparison <- raster('results/10km_esh/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh10km_10perc.tif')
ras_ed_comparison <- raster('results/10km_esh/minshort_speciestargets_waPlants_phylo-ed-comparison__esh10km_10perc.tif')

stopifnot(
  assert_that( compareRaster(ras_threat_comparison,ras_ed_comparison ),compareRaster(ras_threat_comparison,ras_normal_comparison ),
             compareRaster(ras_normal_comparison,land))
)

# Assign a distance class to each in a venn diagramm
ras <- stack(ras_normal_comparison, ras_threat_comparison,ras_ed_comparison )
ras[is.na(ras)] <- 0
names(ras) <- c('normal','threat','ed')
ras <- raster::mask(ras,land) # Mask out land

# Now assign a colour based on whether a certain value is equal/smaller than another
df <- as.data.frame(ras)

# Assign by taxonomic group
df$code <- NA
df$code[which(df$normal>0 & df$threat >0 & df$ed>0)] <- 1 # All layers
df$code[which(df$normal>0 & df$threat == 0 & df$ed==0)] <- 2 # Normal only
df$code[which(df$normal==0 & df$threat >0 & df$ed==0)] <- 3 # Threat only
df$code[which(df$normal==0 & df$threat == 0 & df$ed>0)] <- 4 # ED only
df$code[which(df$normal>0 & df$threat >0 & df$ed==0)] <- 5 # Normal and Threat
df$code[which(df$normal>0 & df$threat ==0 & df$ed>0)] <- 6 # Normal and ED
df$code[which(df$normal==0 & df$threat >0 & df$ed>0)] <- 7 # Threat and ED

# Make output
o <- emptyraster(ras)
o[] <- df$code

# Red = Normal, Green = 
cols = data.frame( code = seq(0,7),
                   hex = c("grey20","#b4d3a6","#ff7f7f","#ffecb2","#87a9e5",
                           "#b49366","#dad686","#6293b9")  
)
# Mask with land again
o[is.na(o)] <- 0
o <- raster::mask(o,land) # Mask out land
        # Save output
writeGeoTiff(o,paste0("results/",target_resolution,"_speciesweights_overlapcode.tif"))

# Make a plot
gm <- ggR(o, layer = 1, maxpixels = 1e10,
         geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  cowplot::theme_map(font_size =  12,font_family = 'Arial') +
  scale_fill_dimanual(values = cols$hex,na.value = 'transparent') + guides(fill = 'none')+
  theme(legend.position = c(.05,.3),legend.background = element_rect(fill = 'transparent')) +
  labs(x = "", y = "" )
gm
ggsave(paste0(figure_path,'/','SIFigure_BiodiversityWeights_map.png'),plot = gm,width = 8,height = 8,dpi = 400)

# How many species in each ?
what = "50km"

phylo_data <- data.frame()
dd <- list.dirs(paste0(results_path,"/",what,target_range,"/"))
if(repr_id) { dd <- dd[grep(repr_run_path,dd)] } else { dd <- dd[1]}
# Get all normal files
ll <- list.files(path = dd,
                 pattern = paste0(
                   'minshort_speciestargets',
                   split_id,
                   '', # PA
                   '', # Carbon
                   '', # Water
                   '', # Carbon Weight
                   '', # Water weight
                   '_waPlants', 
                   '_phylo-noweight-comparison_',
                   target_range,
                   what
                 ),
                 full.names = TRUE
)
ll <- ll[has_extension(ll,'fst')];stopifnot(length(ll)>0)

# All representation budgets
rr <- do.call(rbind, pblapply(ll, function(x,nr = 7) {
  read.fst(x) %>% calc_shortfall %>%
    dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
})
)

# Join in species ids
rr <- rr %>% left_join(.,
                       read_fst(feature_path) %>% mutate(feature = paste0(iucn_id_no,'__',binomial) ) %>%
                         dplyr::select(feature,category, data, kingdom, class) %>% 
                         dplyr::bind_rows( tibble()  ), # Add Carbon!
                       by ='feature' )
# Reformat
rr$budget <- gsub('\\D','',rr$budget)
rr$budget <- factor(rr$budget, levels =seq(10,100,10) )

phylo_data <- bind_rows(phylo_data,
                        rr %>% dplyr::group_by(run,budget) %>% 
  dplyr::summarise(target_reached = sum(target_reached) / n() ) %>% ungroup() %>% 
  dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
  # Average + sd per run
  dplyr::group_by(budget) %>% 
  dplyr::summarise(target_reached_avg = mean(target_reached*100),
                   target_reached_ymin = mean(target_reached*100) - sd(target_reached*100),
                   target_reached_ymax = mean(target_reached*100) + sd(target_reached*100)) %>% 
  dplyr::mutate(type = 'No weight',nr_species = n_distinct(rr$feature))
)

# Weighted by Threat
ll <- list.files(path = dd,
                 pattern = paste0(
                   'minshort_speciestargets',
                   split_id,
                   '', # PA
                   '', # Carbon
                   '', # Water
                   '', # Carbon Weight
                   '', # Water weight
                   '', 
                   '_phylo-rank_weight-comparison_',
                   target_range,
                   what
                 ),
                 full.names = TRUE
)
ll <- ll[has_extension(ll,'fst')];stopifnot(length(ll)>0)

# All representation budgets
rr <- do.call(rbind, pblapply(ll, function(x,nr = 7) {
  read.fst(x) %>% calc_shortfall %>%
    dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
})
)

# Join in species ids
rr <- rr %>% left_join(.,
                       read_fst(feature_path) %>% mutate(feature = paste0(iucn_id_no,'__',binomial) ) %>%
                         dplyr::select(feature,category, data, kingdom, class) %>% 
                         dplyr::bind_rows( tibble()  ), # Add Carbon!
                       by ='feature' )
# Reformat
rr$budget <- gsub('\\D','',rr$budget)
rr$budget <- factor(rr$budget, levels =seq(10,100,10) )

phylo_data <- bind_rows(phylo_data,
                        rr %>% dplyr::group_by(run,budget) %>% 
                          dplyr::summarise(target_reached = sum(target_reached) / n() ) %>% ungroup() %>% 
                          dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
                          # Average + sd per run
                          dplyr::group_by(budget) %>% 
                          dplyr::summarise(target_reached_avg = mean(target_reached*100),
                                           target_reached_ymin = mean(target_reached*100) - sd(target_reached*100),
                                           target_reached_ymax = mean(target_reached*100) + sd(target_reached*100)) %>% 
                          dplyr::mutate(type = 'Threat weight',nr_species = n_distinct(rr$feature))
)

# Weighted by ED
ll <- list.files(path = dd,
                 pattern = paste0(
                   'minshort_speciestargets',
                   split_id,
                   '', # PA
                   '', # Carbon
                   '', # Water
                   '', # Carbon Weight
                   '', # Water weight
                   '_waPlants', 
                   '_phylo-ed-comparison_',
                   target_range,
                   what
                 ),
                 full.names = TRUE
)
ll <- ll[has_extension(ll,'fst')];stopifnot(length(ll)>0)

# All representation budgets
rr <- do.call(rbind, pblapply(ll, function(x,nr = 7) {
  read.fst(x) %>% calc_shortfall %>%
    dplyr::mutate(run = basename(dirname(x)), budget = str_split(tools::file_path_sans_ext( basename(x) ),"_",simplify = TRUE)[,nr] )
})
)

# Join in species ids
rr <- rr %>% left_join(.,
                       read_fst(feature_path) %>% mutate(feature = paste0(iucn_id_no,'__',binomial) ) %>%
                         dplyr::select(feature,category, data, kingdom, class) %>% 
                         dplyr::bind_rows( tibble()  ), # Add Carbon!
                       by ='feature' )
# Reformat
rr$budget <- gsub('\\D','',rr$budget)
rr$budget <- factor(rr$budget, levels =seq(10,100,10) )

phylo_data <- bind_rows(phylo_data,
                        rr %>% dplyr::group_by(run,budget) %>% 
                          dplyr::summarise(target_reached = sum(target_reached) / n() ) %>% ungroup() %>% 
                          dplyr::mutate(budget = as.numeric( as.character( budget ) )) %>% 
                          # Average + sd per run
                          dplyr::group_by(budget) %>% 
                          dplyr::summarise(target_reached_avg = mean(target_reached*100),
                                           target_reached_ymin = mean(target_reached*100) - sd(target_reached*100),
                                           target_reached_ymax = mean(target_reached*100) + sd(target_reached*100)) %>% 
                          dplyr::mutate(type = 'ED weight',nr_species = n_distinct(rr$feature))
)

library(ggforce)
# Construct Plot
phylo_data$type <- factor(phylo_data$type, levels = c('No weight',"Threat weight","ED weight"))

g <- ggplot(phylo_data, aes(x = budget, y = target_reached_avg, group = type, colour = type)) + 
  #theme_few(base_size = 20) +
  theme_bw(base_size = 20) +
  #geom_hline(yintercept = 100,linetype = "dotted", size = .5) +
  geom_ribbon(aes(ymin = target_reached_ymin, ymax = target_reached_ymax), fill = 'grey80', colour =NA, alpha = .7) +
  geom_line(size = 2) +
  scale_color_manual(values = c('red','green','blue'),guide = guide_legend("",override.aes = aes(fill = NA)) ) +
  #scale_linetype_manual(values = c('solid','dotted','twodash'),guide = guide_legend("",override.aes = aes(fill = NA))) + 
  theme(legend.position = c(.25,.15),legend.key.width = unit(.5,'in'),legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent') ) +
  #scale_x_continuous(breaks = pretty_breaks(5),limits = c(as.numeric(levels(rr$budget)[1]),100)) +
  #scale_y_continuous(breaks = pretty_breaks(5),limits = c(0,100)) +
  # Facet zoom:
  facet_zoom(xlim = c(5,15), ylim = c(50, 65),show.area = TRUE,
             horizontal = FALSE) +   theme(zoom.y = element_blank(), validate = FALSE) +
  guides(fill = 'none') +
  ylab("Targets reached (% of all species)") + xlab("Area managed for conservation\n(% of land area)") +
  theme(plot.background = element_blank())
g
ggsave(paste0(figure_path,'/','SIFigure_BiodiversityWeights_plot.png'),plot = g,width = 8,height = 8,dpi = 400)

# ------------------------------------------------------- #
#### [SI-Figure] Plot with range size distributions and plots ####
library(fst)

# Load rij_table
rij_data <- read_fst('data/table_puspecies_50km_esh.fst') %>% dplyr::mutate(amount = as.numeric(amount > 0))
feature_data <- read_fst('data/speciesID_table_50km.fst')

# And the representative sets
lf <- list.files('sets','splist',full.names = T)
repr_id_file <- lapply(lf, function(x) { read_csv(x) %>% dplyr::select(2) } )
names(repr_id_file) <- paste0('id',1:length(repr_id_file))
# Rather quick and dirty rename the columns too
for(i in names(repr_id_file)){
  x = repr_id_file[[i]]; names(x) <- i
  repr_id_file[[i]] <- x
}
rm(lf)

# First recreate the range size distributions
range_sizes <- rij_data %>% dplyr::group_by(id) %>% 
  dplyr::summarise(area = sum(amount,na.rm = TRUE))

# Now for the raster plots
animal_data = feature_data %>%  dplyr::filter(kingdom == 'ANIMALIA')
animal_rij <- rij_data %>% dplyr::filter(id %in% animal_data$iucn_id_no)

# Piero script to create a random plot
rand_animal = data.frame(nrow=3339)
availableids = dplyr::select(animal_data,iucn_id_no)
for (i in 1:10){
  print(i)
  rand_animal=cbind(rand_animal,sample(availableids[,1],size=3339)) # 3339 being approximately 10% of the full dataset
  availableids= as.data.frame(setdiff(availableids[,1], rand_animal[,i+1]))
}
rand_animal=rand_animal[,2:11]
names(rand_animal)=paste('id',1:10,sep='')

df_ranges <- bind_rows(
  range_sizes %>% dplyr::mutate(type = 'full'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[1]]$id1) %>% dplyr::mutate(type = 'set1'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[2]]$id2) %>% dplyr::mutate(type = 'set2'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[3]]$id3) %>% dplyr::mutate(type = 'set3'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[4]]$id4) %>% dplyr::mutate(type = 'set4'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[5]]$id5) %>% dplyr::mutate(type = 'set5'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[6]]$id6) %>% dplyr::mutate(type = 'set6'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[7]]$id7) %>% dplyr::mutate(type = 'set7'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[8]]$id8) %>% dplyr::mutate(type = 'set8'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[9]]$id9) %>% dplyr::mutate(type = 'set9'),
  range_sizes %>% dplyr::filter(id %in% repr_id_file[[10]]$id10) %>% dplyr::mutate(type = 'set10')
)
df_ranges$type <- factor(df_ranges$type,levels = c('full','set1','set2','set3','set4','set5','set6','set7','set8','set9','set10'))  

g_rangesize <- ggplot(df_ranges, aes(x = type, y = log10(area), fill = type)) +
  theme_bw(base_size = 20,base_family = 'Arial') +
  geom_violin(scale = 'area') +
  scale_fill_stata() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = expression(log10('AOH size') ) ) +
  guides(fill = "none")
g_rangesize

# Some statistical test
library(multcomp)
fit <- aov(area ~ type, data = df_ranges %>% mutate(area = log10(area+1)))
summary( glht(fit, linfct = mcp(type = "Tukey") ) )

df_ranges %>% group_by(type) %>%
  summarise(me = median(log10(area+1)),
            m = mean(log10(area+1)))

# plot raster results
require(rgdal)
require(raster)
pu_id_raster  <- raster("data/puid_50km.tif")
land <- raster(paste0("data/globalgrid_mollweide_50km.tif"))
land <- setMinMax(land)
land.mask <- land
land.mask[land.mask>0] <- 1;land.mask[land.mask!=1] <- NA

rich_animal_samp1 = animal_rij %>% dplyr::filter(id %in% repr_id_file$id1$id1) %>% 
  group_by(pu) %>% dplyr::summarise(sr = n_distinct(id) ) %>% 
  dplyr::select( pu, sr ) %>% 
  drop_na()

richness<-function(rijdata,colnum){
  rijdata %>% dplyr::filter(id %in% rand_animal[,colnum]) %>% 
    group_by(pu) %>% dplyr::summarise(sr = n_distinct(id) ) %>% 
    dplyr::select( pu, sr ) %>% 
    drop_na()
}

richness2<-function(rijdata,colnum){
  rijdata %>% dplyr::filter(id %in% repr_id_file[[paste0("id",colnum)]][[1]] ) %>% 
    group_by(pu) %>% dplyr::summarise(sr = n_distinct(id) ) %>% 
    dplyr::select( pu, sr ) %>% 
    drop_na()
}

rich_animal_samp2<-richness(animal_rij,2)
rich_animal_samp3<-richness(animal_rij,3)
rich_animal_samp4<-richness(animal_rij,4)
rich_animal_samp5<-richness(animal_rij,5)
rich_animal_samp6<-richness(animal_rij,6)
rich_animal_samp7<-richness(animal_rij,7)
rich_animal_samp8<-richness(animal_rij,8)
rich_animal_samp9<-richness(animal_rij,9)
rich_animal_samp10<-richness(animal_rij,10)

rich_animal_powo1<-richness2(animal_rij,1)
rich_animal_powo2<-richness2(animal_rij,2)
rich_animal_powo3<-richness2(animal_rij,3)
rich_animal_powo4<-richness2(animal_rij,4)
rich_animal_powo5<-richness2(animal_rij,5)

  
rich_animal_total=animal_rij  %>% 
  group_by(pu) %>% dplyr::summarise(sr = n_distinct(id) ) %>% 
  dplyr::select( pu, sr ) %>% 
  drop_na()

animal_rich_map <- createOutput(pu_id_raster, as.data.frame(rich_animal_total), scol = 'sr' )
animal_rich_map_sub2 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp3), scol = 'sr' )
plot(animal_rich_map_sub2)
animal_rich_map_sub3 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp3), scol = 'sr' )
animal_rich_map_sub4 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp4), scol = 'sr' )
animal_rich_map_sub5 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp5), scol = 'sr' )
animal_rich_map_sub6 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp6), scol = 'sr' )
animal_rich_map_sub7 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp7), scol = 'sr' )
animal_rich_map_sub8 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp8), scol = 'sr' )
animal_rich_map_sub9 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp9), scol = 'sr' )
animal_rich_map_sub10 <- createOutput(pu_id_raster, as.data.frame(rich_animal_samp10), scol = 'sr' )

animal_rich_map_powo1 <- createOutput(pu_id_raster, as.data.frame(rich_animal_powo1), scol = 'sr' )
plot(animal_rich_map_powo1)
animal_rich_map_powo2 <- createOutput(pu_id_raster, as.data.frame(rich_animal_powo2), scol = 'sr' )
animal_rich_map_powo3 <- createOutput(pu_id_raster, as.data.frame(rich_animal_powo3), scol = 'sr' )
animal_rich_map_powo4 <- createOutput(pu_id_raster, as.data.frame(rich_animal_powo4), scol = 'sr' )
animal_rich_map_powo5 <- createOutput(pu_id_raster, as.data.frame(rich_animal_powo5), scol = 'sr' )

# Normalize them
results <- raster::stack()
for(i in 2:10){
  results <- addLayer(results, normalize_raster(get(paste0('animal_rich_map_sub',i))) )
}
results_powo <- raster::stack()
for(i in 1:5){
  results_powo <- addLayer(results_powo, normalize_raster(get(paste0('animal_rich_map_powo',i))) )
}

norm_animal_rich_map <- normalize_raster(animal_rich_map)
plot(norm_animal_rich_map - mean(results_powo))

pr <- stack(mean(results),mean(results_powo),norm_animal_rich_map )
names(pr) <- c("Random","wgsrpd","Full")
plot(pr)
pr[is.na(pr)] <- 0
pr <- raster::mask(pr, land.mask )

# Do a statistical test
ex <- as.data.frame(stack(pr,pu_id_raster)) %>% drop_na()
cor.test(ex$Full,ex$puid_50km,method = "kendall")
cor.test(ex$Full,ex$Random,method = "kendall")
cor.test(ex$Full,ex$wgsrpd,method = "kendall")

library(boot)
set.seed(1)
bb <- boot(ex, 
           statistic = function(data, i) {
             cor(data[i, "Full"], data[i, "Random"], method='pearson')
           },
           R = 1000
)
bb

library(RStoolbox)
names(pr)
g_maps <- ggR(pr,layer = c(3,1,2),geom_raster = TRUE) +
  ggthemes::theme_map(base_size = 20,base_family = 'Arial') +
#  theme_bw(base_size = 20,base_family = 'Arial') +
  scale_fill_gradientn(colours = scico(10,palette = 'roma',direction = -1),na.value = 'transparent') +
  guides(fill = 'none')
g_maps

library(cowplot)
pg <- plot_grid(g_maps,g_rangesize,nrow = 2,labels = 'auto')
cowplot::ggsave2(filename = 'figures/SIFigure_SamplingRandom.png',width = 12,height = 8,dpi = 400)

#### [SI-Table] National ranking for high priority regions ####
# National responsibility for global areas of significance
# Also with currently protected areas
library(raster)
library(sf)
library(units)
library(exactextractr)
library(fasterize)
rasterOptions(progress = 'text')

# Create a national ranking of priority regions for each country
# -> Extract values from the biodiv only, biodiv + carbon, biodiv + carbon + water from the map
# -> Calculate an area-weighted score from these values (e.g. multiplied by the proportional abundance of terrestrial land area [i.e., area (km2) divided by total area] )
# -> Export table

# Load priority ranked raster
ras_biodivonly <- raster(paste0(results_path,"/","10km_esh_ranked/minshort_speciestargets_esh10km_repruns10_ranked.tif"))
ras_biocarbon <- raster(paste0(results_path,"/","10km_esh_ranked/minshort_speciestargets_carbon__esh10km_repruns10_ranked.tif"))
ras_biowater <- raster(paste0(results_path,"/","10km_esh_ranked/minshort_speciestargets_water__esh10km_repruns10_ranked.tif"))
ras_biocarbonwater <- raster(paste0(results_path,"/","10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked.tif"))
ras_land <- raster('data/globalgrid_mollweide_10km.tif')
stopifnot(assert_that(compareRaster(ras_land,ras_biodivonly)))

# For countries use Natural Earth
countries <- st_read('/media/martin/data/raw/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp') %>% 
  dplyr::select(NAME,ABBREV,ECONOMY,INCOME_GRP,CONTINENT,SUBREGION,REGION_WB) %>% 
  dplyr::filter(CONTINENT %notin% 'Antarctica') %>% 
  st_transform(crs(proj4string(ras_biodivonly))) %>% # Transform to mollweide
  # Calculate area
  dplyr::mutate(area = sf::st_area(.) ) %>% 
  # Proportional area
  dplyr::mutate(prop_area = (area /  sum(area)) )

df.countries <- as.data.frame(countries %>% st_drop_geometry() %>% dplyr::select(NAME))

# Fasterize the proportional area as weights
weights_raster <- fasterize(countries, raster = ras_land,field = "prop_area",fun = 'max')

# Now extract the weighted priority estimate per country globally. Multiply with proprotional 
#ex <- raster::extract(ras_priority, countries, fun = sum, na.rm= TRUE, small = TRUE, weights = TRUE,normalizeWeights = TRUE, df = TRUE)
df.countries$m_biodivonly <- exactextractr::exact_extract(x = ras_biodivonly, y = countries, fun = "mean", progress = TRUE)
df.countries$wm_biodivonly <- exactextractr::exact_extract(x = ras_biodivonly*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.countries$m_biocarbon <- exactextractr::exact_extract(x = ras_biocarbon, y = countries, fun = "mean", progress = TRUE)
df.countries$wm_biocarbon <- exactextractr::exact_extract(x = ras_biocarbon*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.countries$m_biowater <- exactextractr::exact_extract(x = ras_biowater, y = countries, fun = "mean", progress = TRUE)
df.countries$wm_biowater <- exactextractr::exact_extract(x = ras_biowater*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.countries$m_biocarbonwater <- exactextractr::exact_extract(x = ras_biocarbonwater, y = countries, fun = "mean", progress = TRUE)
df.countries$wm_biocarbonwater <- exactextractr::exact_extract(x = ras_biocarbonwater*weights_raster, y = countries, fun = "mean", progress = TRUE)
# Drop those that got no values (too small, islands, etc...)
df.countries <- df.countries %>% drop_na(m_biocarbonwater)
# Finally rank them all
df.countries2 <- df.countries %>% mutate_at(c("m_biodivonly","m_biocarbon","m_biowater","m_biocarbonwater",
                                             "wm_biodivonly","wm_biocarbon","wm_biowater","wm_biocarbonwater"), rank) %>% 
  # Sort
  arrange(m_biodivonly)

# Rename for export
export1 <- df.countries2 %>% rename(
  "Country name" = "NAME",
  "Biodiversity rank (mean)" = "m_biodivonly",
  "Biodiversity_Carbon rank (mean)" = "m_biocarbon",
  "Biodiversity_Water rank (mean)" = "m_biowater",
  "Biodiversity_Carbon_Water rank (mean)" = "m_biocarbonwater",
  "Biodiversity rank (area-weighted mean)" = "wm_biodivonly",
  "Biodiversity_Carbon rank (area-weighted mean)" = "wm_biocarbon",
  "Biodiversity_Water rank (area-weighted mean)" = "wm_biowater",
  "Biodiversity_Carbon_Water rank (area-weighted mean)" = "wm_biocarbonwater"
)
write_excel_csv2(x = export1, "SITable3_ImportanceRanking_country.csv",row.names = FALSE)

# Load biomes shapefiles
biomes_code <- sf::read_sf('/media/martin/data/raw/Ecoregions2017/Ecoregions2017.shp') %>% 
  dplyr::select(BIOME_NUM,BIOME_NAME,COLOR_BIO) %>% 
  dplyr::arrange(BIOME_NUM) %>% group_by(BIOME_NUM,BIOME_NAME) %>% summarise()

# Also do per biome
df.biome$m_biodivonly <- exactextractr::exact_extract(x = ras_biodivonly, y = countries, fun = "mean", progress = TRUE)
df.biome$wm_biodivonly <- exactextractr::exact_extract(x = ras_biodivonly*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.biome$m_biocarbon <- exactextractr::exact_extract(x = ras_biocarbon, y = countries, fun = "mean", progress = TRUE)
df.biome$wm_biocarbon <- exactextractr::exact_extract(x = ras_biocarbon*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.biome$m_biowater <- exactextractr::exact_extract(x = ras_biowater, y = countries, fun = "mean", progress = TRUE)
df.biome$wm_biowater <- exactextractr::exact_extract(x = ras_biowater*weights_raster, y = countries, fun = "mean", progress = TRUE)
df.biome$m_biocarbonwater <- exactextractr::exact_extract(x = ras_biocarbonwater, y = countries, fun = "mean", progress = TRUE)
df.biome$wm_biocarbonwater <- exactextractr::exact_extract(x = ras_biocarbonwater*weights_raster, y = countries, fun = "mean", progress = TRUE)

# ------------------------------------------------------- #
#### [Response figure] Ranking of individual sets ####
# Script to create ranks of the representative sets individually

assert_that(exists('ll.tif'))
newdir <- 'results/10km_esh_ranked_sets'
dir.create(newdir,showWarnings = FALSE)

plot.list <- list()
# Loop through each individual set and create a ranked version
for(set in paste0('id',1:10) ){
  print(set)
  # Get the individual set
  ll.tif.set <- ll.tif[grep(set,ll.tif)]
  if(length(ll.tif.set)>10) ll.tif.set <- ll.tif.set[1:10] # For lvl 1 that get matched with 10
  
  # Load raster
  ras2 <- makeStack(ll.tif.set)
  #tifs <- grep(c('50perc'),ll.tif,invert = F,value = T)
  ras2 <- reclassify(ras2, cbind(NA, NA, 0), right = FALSE)
  
  ## Mean stack
  mcp_ilp_hier <- mean(ras2, na.rm = T) # Mean! 
  # Mask with global land area mask
  mcp_ilp_hier <- raster::mask(mcp_ilp_hier, land.mask)
  
  mcp_ilp_hier_norm <- mcp_ilp_hier * 100 # precision of input layers
  # -- #
  # Remove 0 fractional estimates
  mcp_ilp_hier_norm[mcp_ilp_hier_norm==0] <- NA
  # Convert to ranks and invert
  mcp_ilp_hier_norm <- abs( raster_to_ranks(mcp_ilp_hier_norm,n = 100,plot = FALSE) - 100 )
  # Mask again with global land area mask 
  mcp_ilp_hier_norm[is.na(mcp_ilp_hier_norm)] <- 100 # Add lest important areas as 100 so that proportions match again
  mcp_ilp_hier_norm <-  raster::mask(mcp_ilp_hier_norm,land.mask)
  NAvalue(mcp_ilp_hier_norm) <- -9999
  
  (fname <- paste0(newdir,"/minshort_speciestargets",
                   split_id,pa_fname,pa_carbon,pa_water,pa_carbmult,pa_watermult,pa_plants,target_range,target_resolution,
                   '_',set,
                   "_ranked.tif"))
  # Save output
  writeGeoTiff(mcp_ilp_hier_norm,fname)
  
  g <- ggR(mcp_ilp_hier_norm, layer = 1, maxpixels = 1e10,
           geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
    ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
    scale_fill_gradientn(breaks = c(100,75,50,25,1),colours = scico(10,palette = 'roma'),na.value = 'transparent',
                         guide = guide_colourbar(title = "Priority rank",title.theme = element_text(size = 14),
                                                 title.position = 'bottom',title.hjust = 0.5,label.position = 'bottom',barwidth = unit(1.5,'in'),
                                                 direction = 'horizontal',ticks = F,reverse = T) ) +
    theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent')) +
    guides(fill = 'none',colour = 'none') +
    labs(x = "", y = "", title = paste0('Set: ', set) ) +
    theme(title = element_text(size = 16),plot.title = element_text(hjust = .5)) +
    theme(plot.margin = unit(c(0,0,0,0),'cm'))
  ggsave(filename = paste0(tools::file_path_sans_ext(fname), '.png'),
         plot = g,width = 10,height = 6,units = 'in',dpi = 400)
  
  # Save in list
  plot.list[[set]] <- g
}

# Now get all plots together
library(cowplot)
library(gridExtra)
pg <- cowplot::plot_grid(plotlist = plot.list,
                         nrow = 3,ncol = 4)
cowplot::ggsave2(plot = pg,filename = 'figures/ResponseFigure_RankedSets.png',width = 14,height = 14,units = 'in')
