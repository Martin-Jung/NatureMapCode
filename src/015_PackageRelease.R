library(raster)
library(stringr)
source('src/000_ConvenienceFunctions.R')

createFolder <- function(name) dir.create(name,showWarnings = FALSE)
makePlot <- function(name,fname, cols = scico::scico(10,palette = 'roma') ){
  r <- raster(name)
  fname <- paste0(tools::file_path_sans_ext(fname),'.png')
  if(stringr::str_detect(name, 'cv')) cols = scico(10,palette = 'tokyo',direction = -1)
  png(fname, width = 1000, height = 800, units = 'px', type = 'cairo')
  plot(r, col = cols)
  dev.off()
}

# Define output path
# path_output = 'P:/bec_data/001_projects/2019_NatureMap'
path_output = '/mnt/pdrive/bec_data/001_projects/2019_NatureMap'
createFolder(path_output)

# Define output folder
path_output_folder = paste0(path_output,'/NatureMap_prioritymaps')
createFolder(path_output_folder)

# ------------- #
##############################
# Save the global grid ids and fractional land cover data
land <- raster("data/globalgrid_mollweide_10km.tif")
land <- setMinMax(land)
writeGeoTiff(land, paste0(path_output_folder,'/globalgrid_mollweide_10km.tif'))
# And 50km
land <- raster("data/globalgrid_mollweide_50km.tif")
land <- setMinMax(land)
writeGeoTiff(land, paste0(path_output_folder,'/globalgrid_mollweide_50km.tif'))
##############################

#### 'BiodiversityCarbonWater' ----
newdir = paste0(path_output_folder,'/BiodiversityCarbonWater/10km')
createFolder(newdir)
# Copy files over
file.copy(from = "results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_ranked.tif",
          to = paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_ranked.tif")) # Mean
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_ranked.tif",
         paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_ranked.png"))
file.copy(from = "results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_cv.tif",
          to = paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_cv.tif")) # Mean
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_cv.tif",
         paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh10km_repruns10_cv.png"),
         cols = scico(10,palette = 'tokyo',direction = -1))

# --- #
# With Protected areas locked in
file.copy('results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_ranked.tif',
          paste0(newdir, '/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_ranked.tif'))
makePlot('results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_ranked.tif',
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_ranked.png'))

file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_cv.tif",
          paste0(newdir, '/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_cv.tif'))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__water__esh10km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# --------------------- #
# Without biome targets
file.copy(from = "results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked.tif",
          to = paste0(newdir, "/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked.tif")) # Mean
makePlot("results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked.tif",
         paste0(newdir, "/minshort_speciestargets_carbon__water__esh10km_repruns10_ranked.png"))
file.copy(from = "results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_cv.tif",
          to = paste0(newdir, "/minshort_speciestargets_carbon__water__esh10km_repruns10_cv.tif")) # Mean
makePlot("results/10km_esh_ranked/minshort_speciestargets_carbon__water__esh10km_repruns10_cv.tif",
         paste0(newdir, "/minshort_speciestargets_carbon__water__esh10km_repruns10_cv.png"),
         cols = scico(10,palette = 'tokyo',direction = -1))
# 
# # --- #
# # Without Biome targets and Protected areas locked in
file.copy('results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh10km_repruns10_ranked.tif',
          paste0(newdir, '/minshort_speciestargetswithPA_carbon__water__esh10km_repruns10_ranked.tif'))
makePlot('results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh10km_repruns10_ranked.tif',
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__water__esh10km_repruns10_ranked.png'))

# --- # 
# And for 50km 
newdir = paste0(path_output_folder,'/BiodiversityCarbonWater/50km')
createFolder(newdir)

file.copy("results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_carbon__water__esh50km_repruns10_ranked.png'))

file.copy("results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargets_carbon__water__esh50km_repruns10_cv.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_carbon__water__esh50km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_carbon__water__esh50km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# With PA
file.copy("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_ranked.png'))

file.copy("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_cv.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__water__esh50km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# With Biome id and PA
file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_ranked.png'))

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_cv.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__water__esh50km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# Without PA and biome id
file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_ranked.png'))

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_cv.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__carbon__water__esh50km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))


# ---------------- #
#### 'BiodiversityCarbon' ----
createFolder(paste0(path_output_folder,'/BiodiversityCarbon/'))
newdir = paste0(path_output_folder,'/BiodiversityCarbon/10km')
createFolder(newdir)

# Without Biome id and without PA
file.copy("results/10km_esh_ranked/minshort_speciestargets_carbon__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_carbon__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_carbon__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_carbon__esh10km_repruns10_ranked.png'))

# Without Biome id and with PA
file.copy("results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_carbon__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__esh10km_repruns10_ranked.png'))

file.copy("results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__esh10km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_carbon__esh10km_repruns10_cv.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargetswithPA_carbon__esh10km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__esh10km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# With Biome ID and without PA
file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_ranked.png'))

file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_cv.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__carbon__esh10km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# With Biome ID and with PA
file.copy( "results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_ranked.tif"))
makePlot( "results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_ranked.png'))

file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_cv.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_cv.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_cv.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__esh10km_repruns10_cv.png'),
         cols = scico(10,palette = 'tokyo',direction = -1))

# --- #
# 50km
newdir = paste0(path_output_folder,'/BiodiversityCarbon/50km')
createFolder(newdir)
file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_carbon__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_carbon__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_carbon__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__carbon__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__carbon__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id__carbon__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_carbon__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_carbon__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_carbon__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_carbon__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_carbon__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_carbon__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_carbon__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargetswithPA_carbon__esh50km_repruns10_cv.tif"

# ---------------- #
#### 'BiodiversityWater' ----
createFolder(paste0(path_output_folder,'/BiodiversityWater'))
newdir = paste0(path_output_folder,'/BiodiversityWater/10km')
createFolder(newdir)

# No PA and without biome
file.copy("results/10km_esh_ranked/minshort_speciestargets_water__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_water__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_water__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_water__esh10km_repruns10_ranked.png'))
# "results/10km_esh_ranked/minshort_speciestargets_water__esh10km_repruns10_cv.tif"

# PA and without biome
file.copy("results/10km_esh_ranked/minshort_speciestargetswithPA_water__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_water__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargetswithPA_water__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_water__esh10km_repruns10_ranked.png'))
# "results/10km_esh_ranked/minshort_speciestargetswithPA_water__esh10km_repruns10_cv.tif"

# PA and with biome
file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_water__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_water__esh10km_repruns10_ranked.png'))
# "results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh10km_repruns10_cv.tif"

# NO Pa and with biome
file.copy( "results/10km_esh_ranked/minshort_speciestargets_biome.id__water__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__water__esh10km_repruns10_ranked.tif"))
makePlot( "results/10km_esh_ranked/minshort_speciestargets_biome.id__water__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__water__esh10km_repruns10_ranked.png'))
# "results/10km_esh_ranked/minshort_speciestargets_biome.id__water__esh10km_repruns10_cv.tif"

# --- #
# 50km
newdir = paste0(path_output_folder,'/BiodiversityWater/50km')
createFolder(newdir)

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id__water__esh50km_repruns10_ranked.tif",
           paste0(newdir, "/minshort_speciestargets_biome.id__water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id__water__esh50km_repruns10_ranked.tif",
          paste0(newdir,'/minshort_speciestargets_biome.id__water__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id__water__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_water__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_water__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_water__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_water__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargetswithPA_water__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_water__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargetswithPA_water__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_water__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargetswithPA_water__esh50km_repruns10_cv.tif"

# ---------------- #
#### 'Biodiversityonly' ----
createFolder(paste0(path_output_folder,'/BiodiversityOnly'))
newdir = paste0(path_output_folder,'/BiodiversityOnly/10km')
createFolder(newdir)

file.copy("results/10km_esh_ranked/minshort_speciestargets_esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_esh10km_repruns10_ranked.png'))
#"results/10km_esh_ranked/minshort_speciestargets_esh10km_repruns10_cv.tif"

# Only 1 set
file.copy("results/10km_esh_ranked/minshort_speciestargetswithPA_esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargetswithPA_esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_esh10km_repruns10_ranked.png'))
#"results/10km_esh_ranked/minshort_speciestargetswithPA_esh10km_repruns10_cv.tif"

file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_esh10km_repruns10_ranked.png'))
#"results/10km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh10km_repruns10_cv.tif"

file.copy("results/10km_esh_ranked/minshort_speciestargets_biome.id__esh10km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__esh10km_repruns10_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_biome.id__esh10km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__esh10km_repruns10_ranked.png'))
# "results/10km_esh_ranked/minshort_speciestargets_biome.id__esh10km_repruns10_cv.tif"

# ------- #
# 50KM
newdir = paste0(path_output_folder,'/BiodiversityOnly/50km')
createFolder(newdir)

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id__esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id__esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id__esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id__esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id__esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargets_esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargets_esh50km_repruns10_cv.tif"

file.copy("results/50km_esh_ranked/minshort_speciestargetswithPA_esh50km_repruns10_ranked.tif",
          paste0(newdir, "/minshort_speciestargetswithPA_esh50km_repruns10_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargetswithPA_esh50km_repruns10_ranked.tif",
         paste0(newdir,'/minshort_speciestargetswithPA_esh50km_repruns10_ranked.png'))
#"results/50km_esh_ranked/minshort_speciestargetswithPA_esh50km_repruns10_cv.tif"

# -------------- #
#### Phylogenetic ####
createFolder(paste0(path_output_folder,'/BiodiversityWeights'))
newdir = paste0(path_output_folder,'/BiodiversityWeights/50km')
createFolder(newdir)

# No weights
file.copy("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh50km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh50km_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh50km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh50km_ranked.png'))

# Rank weights
file.copy("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh50km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh50km_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh50km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh50km_ranked.png'))

# ED weighted
file.copy("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-ed-comparison__esh50km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-ed-comparison__esh50km_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-ed-comparison__esh50km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-ed-comparison__esh50km_ranked.png'))

# EDGE weighted
file.copy("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-edge-comparison__esh50km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-edge-comparison__esh50km_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-edge-comparison__esh50km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-edge-comparison__esh50km_ranked.png'))

# --- #
# 50km
createFolder(paste0(path_output_folder,'/BiodiversityWeights'))
newdir = paste0(path_output_folder,'/BiodiversityWeights/50km')
createFolder(newdir)

# No weights
file.copy("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh10km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh10km_ranked.tif"))
makePlot("results/50km_esh_ranked/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh10km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-noweight-comparison__esh10km_ranked.png'))

# Rank weights
file.copy("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh10km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh10km_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh10km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-rank_weight-comparison__esh10km_ranked.png'))

# ED weighted
file.copy("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-ed-comparison__esh10km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-ed-comparison__esh10km_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-ed-comparison__esh10km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-ed-comparison__esh10km_ranked.png'))

# EDGE weighted
file.copy("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-edge-comparison__esh10km_ranked.tif",
          paste0(newdir, "/minshort_speciestargets_waPlants_phylo-edge-comparison__esh10km_ranked.tif"))
makePlot("results/10km_esh_ranked/minshort_speciestargets_waPlants_phylo-edge-comparison__esh10km_ranked.tif",
         paste0(newdir,'/minshort_speciestargets_waPlants_phylo-edge-comparison__esh10km_ranked.png'))
