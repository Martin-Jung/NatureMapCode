# Code for data preparation, analysis and visualization for the NatureMap analysis

This repository contains the code to showcase the steps done in the analysis by Jung et al. and is separated into data preparation scripts, analysis scripts and visualization scripts.  

**Usage Notes**:   
* The purpose of this repository is mainly to show interested third parties the analysis steps undertaken.  
* Scripts contain links to folders which are system-specific (e.g. data storage) and are thus not necessarily reproducible as such.  
* Owing to data-sharing agreement it is not possible to share any raw or intermediate result created in this project. The 'data' folder is thus empty.  
* Created results as part of this project will be uploaded separately to a repository with a digital object identifier (DOI). The 'results' folder is thus empty. See final manuscript for the link to the data respository.  
  
**Preliminary Citation**:  

> Jung et al. (2020) bioRxiv 2020.04.16.021444; doi: https://doi.org/10.1101/2020.04.16.021444 

<hr>

## Purpose of scripts

+ [000_ConvenienceFunctions.R](src/000_ConvenienceFunctions.R) contains supplementary analysis and visualization functions to be loaded at each script.  
+ [001_CreateGlobalFishnets.R](src/001_CreateGlobalFishnets.R) contains functions to create a global background layer of planning units based either GADM or NaturalEarth.  
+ [002_PrepCarbonData.R](src/002_PrepCarbonData.R) contains data prepration functions and code to prepare the carbon data for the prioritization analysis.  
+ [002_PrepWaterData.R](src/002_PrepWaterData.R) contains data prepration functions and code to prepare the water data for the prioritization analysis.  

+ [002_PrepConservationFeatures_BIEN.R](src/002_PrepConservationFeatures_BIEN.R) contains data prepration functions and code to prepare the BIEN species data for the prioritization analysis.  
+ [002_PrepConservationFeatures_BIENAustralia.R](src/002_PrepConservationFeatures_BIENAustralia.R) contains data prepration functions and code to prepare the BIEN Australian species data for the prioritization analysis. + [002_PrepConservationFeatures_GARD.R](src/002_PrepConservationFeatures_GARD.R) contains data prepration functions and code to prepare the GARD species data for the prioritization analysis.  
+ [002_PrepConservationFeatures_IUCN.R](src/002_PrepConservationFeatures_IUCN.R) contains data prepration functions and code to prepare the IUCN species data for the prioritization analysis.  
+ [002_PrepConservationFeatures_KewPlants.R](src/002_PrepConservationFeatures_KewPlants.R) contains data prepration functions and code to prepare the KEW plant species data for the prioritization analysis.  
+ [002_PrepConservationFeatures_NewBIENModels.R](src/002_PrepConservationFeatures_NewBIENModels.R) contains data prepration functions and code to prepare the new BIEN modelled species data for the prioritization analysis.  

+ [003_ArtificalMasking.R](src/003_ArtificalMasking.R) contains functions to create masked AOH from the plant data only using a subset of the global habitat type map.  
+ [003_LoadAndPrepESH.R](src/003_LoadAndPrepESH.R) contains functions to prepare the extracted AOH from Google Earth Engine and place them into individual raster layers.  
+ [003_PrepConstraints.R](src/003_PrepConstraints.R) prepares constrains such as the WDPA dataset and aligns it with the Planning units.  

+ [004_CreateSpeciesMatchupTable.R](src/004_CreateSpeciesMatchupTable.R) creates a full species matchup table to link species names and ids in the prioritization. To be run before the species-pu table is created.    

+ [005_CreateGlobalPUTable.R](src/005_CreateGlobalPUTable.R) Loads all previously prepared layers and creates a global species-pu-amount table (or rij-matrix) containing the amount of each asset per planning unit.  
+ [006_PrepWeights.R](src/006_PrepWeights.R) Formats species-specific weights such as ED or EDGE scores for the analysis.  


+ [008_GlobalProblems.R](src/008_GlobalProblems.R) Main analysis script that creates the global prioritizations dependend on chosen parameters. Returns a file with the solution and shortfall for each budget and file.  
+ [009_CarbonWaterWeightProblems.R](src/009_CarbonWaterWeightProblems.R) Script that is sourced from 008_GlobalProblems to create a prioritization with varying carbon and water weights.  
+ [009_CarbonWeightProblems.R](src/009_CarbonWeightProblems.R) Script that is sourced from 008_GlobalProblems to create a prioritization with varying biodiversity and carbon weights.  


+ [011_FiguresSummary.R](src/011_FiguresSummary.R) Main script to create ranked maps, all figures and summary statistics used in the manuscript based on the files created in 008_GlobalProblems.R.  

+ [015_PackageRelease.R](src/015_PackageRelease.R) Script to prepare all created files for release.  