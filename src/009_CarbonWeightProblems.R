stopifnot( assert_that(exists('rij_data'), exists('feature_data')) ,
           include_carbon == TRUE & include_water == FALSE) # Water not yet implemented!

# Create a new results_path for this problem
results_path <- paste0('results/',target_resolution,target_range,'_',run_variety)
if(!dir.exists(results_path)) { dir.create(results_path) }

# Set a number of weights for each budget target 
target_weight <- round( seq(1, carbonwater_weights[length(carbonwater_weights)],length.out = 10) )

myLog('Loaded carbon_weights problem data for ', n_distinct(rij_data$pu), ' PUs',' and ', n_distinct(rij_data$id), ' features')
# ------------------------------------------------------------------ #
set_number_of_threads( cores ) # Set parallel processing

if(pa_lockedin) {
  budgets <- sort( c(min_budget, budgets[which(budgets > min_budget)] ) ) 
  pu_data$cost <- pu_data$cost / 1000 ; pu_data$protected <- pu_data$protected / 1000
  rij_data$amount <- rij_data$amount / 1000
} else {budgets <- sort( budgets ) }


# Loop through each weight
for(carbon_multiplier in target_weight){
  carbonwater_weights[length(carbonwater_weights)]  <- carbon_multiplier
  
  # and through each budget
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
    pa_carbmult <- paste0("_carbweight",carbon_multiplier,"_")
    # Is there a water weight multiplier greater than 1
    pa_watermult <- ""#ifelse(water_multiplier > 1, paste0("_waterweight",water_multiplier,"_"),"")
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
    # Lock in previous one if higher than 1
    if(k>1){
      assert_that(file.exists(out_previous_save))
      prev_solution <- read_fst(out_previous_save)
      prev_solution <- subset(prev_solution,select = c("id","solution_1")) %>% dplyr::rename(pu = "id",lower = "solution_1") %>% 
        dplyr::filter(lower > 0) %>% # Filter out unselected PUs
        dplyr::mutate(upper = 1) %>% # Add upper bound of whatever the maximum is 
        dplyr::mutate(lower = pmin(lower, 1 )) # Assert that lower does not overshoot 1
      # Lock in previous solution if it was binary or as continious
      stopifnot( assert_that(nrow(prev_solution)>0 ) )
    }
    
    # Budget formulation
    b_cells <- b * sum(pu_data$cost) #n_distinct(pu_data$id)
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
    if(phylo_weights & phylo_comparison == FALSE){ p_maxfeat <- p_maxfeat %>% add_feature_weights(phylo_score) }
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
    out <- feature_representation2(p_maxfeat,r_maxfeat[, "solution_1", drop = FALSE], cores)
    # Join in amount necessary for reaching target as well as feature abundance in planning_units()
    out$absolute_target <- p_maxfeat$targets$output()$value
    out$feature_abundance_pu <- p_maxfeat$feature_abundances_in_planning_units()
    
    # Append target
    write.fst(out, out_name )
    
    # Create raster output
    out_ras <- createOutput(pu_id_raster,r_maxfeat,"solution_1")
    names(out_ras) <- paste0("percent",round(b,2) * 100) # Rename
    # Save output
    writeGeoTiff(out_ras, out_name_tif,dt = ifelse(decision_type=="binary","INT2S","FLT4S") )
    
    # Finally also do this for any existing portfolios / comparable optimal solution
    n_solutions <- grep('solution', names(r_maxfeat),value = TRUE)
    if(port & length(n_solutions) > 1 ){
      myLog('Multiple optimal solutions found. Sampling rasters and representation')
      
      pool_solutions <- data.frame()
      pool_stack <- raster::stack()
      for(sol in n_solutions){
        pool_solutions <- bind_rows(
          pool_solutions,
          feature_representation2(p_maxfeat,r_maxfeat[, sol, drop = FALSE]) %>% 
            dplyr::mutate(
              absolute_target = p_maxfeat$targets$output()$value,
              feature_abundance_pu = p_maxfeat$feature_abundances_in_planning_units(),
              solution = sol
            )
        )
        # Add solution to stack
        pool_stack <- addLayer(pool_stack, createOutput(pu_id_raster,r_maxfeat,sol) )
      }
      # Write pool outputs to new 
      write.fst(out, paste0(results_path,"/","pool_",tools::file_path_sans_ext(x = basename(out_name)),".fst") )
      writeGeoTiff(out_ras, paste0(results_path,"/","pool_",tools::file_path_sans_ext(x = basename(out_name_tif)),".tif"),dt = ifelse(decision_type=="binary","INT2S","FLT4S") )
    }
    rm(r_maxfeat,p_maxfeat,out);gc()
  }
  
}
stop("DONE!")
