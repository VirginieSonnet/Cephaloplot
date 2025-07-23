#' ============================= FUNCTIONS =================================
#' Functions to extract data from model output from Cephalopod. 
#' 
#' * extract_observations: retrieve the original observations 
#' * extract_performance: retrieve predictive performance and cumulative variable importance per model 
#' * extract_predictors: retrieve the input predictors or the predictors selected 
#' * extract_predictions: retrieve the projected observations 
#' =============================================================================




# ----- extract_observations -----

extract_observations <- function(project_wd,FOLDER_NAME,SUBFOLDER_NAME){
  # --- Function to extract the initial observations ---
  # Inputs: 
  #       - project_wd: string, path to the project folder
  #       - FOLDER_NAME: string, folder output from Cephalopod 
  #       - SUBFOLDER_NAME: string, worms_id/scientificname of selected species
  # Outputs: 
  #       - observations: tibble, initial measurementvalue, latitude, longitude, month, year and varfactor (if exist)
  
  
  # 1. Load the query output
  load(paste0(project_wd, FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # 2. Retrieve the measurementvalue
  observations <- bind_cols(QUERY$Y,QUERY$S)
  
  # 3. Return output 
  return(observations)
}


# ----- extract_performance -----

extract_performance <- function(project_wd,FOLDER_NAME,SUBFOLDER_NAME,
                                model=NULL){
  # --- Function to extract the initial observations ---
  # Inputs: 
  #       - project_wd: string, path to the project folder
  #       - FOLDER_NAME: string, folder output from Cephalopod 
  #       - SUBFOLDER_NAME: string, worms_id/scientificname of selected species
  #       - model: vector (optional), specific model(s) to retrieve performance for
  # Outputs: 
  #       - observations: tibble, initial measurementvalue, latitude, longitude, month, year and varfactor (if exist)
  
  
  # 1. Load the model output
  load(paste0(project_wd, FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData")) 
  
  # 2. Compute the models to extract from
  # models that were ran
  available_models <- names(MODEL)[!(names(MODEL) %in% c("MODEL_LIST","recommandations"))]
  
  # if model is NULL, use all available
  models_to_use <- if (is.null(model)) {
    available_models
  } else {
    intersect(available_models, model)
  }
  
  # 3. Retrieve the performance
  performance <- map_dfr(models_to_use, function(mod_name) {
    perf <- as_tibble_row(MODEL[[mod_name]]$eval) %>% 
      mutate(model = mod_name)
  })
  
  # 4. Return output 
  return(performance)
}



# ----- extract_predictors -----

extract_predictors <- function(project_wd,FOLDER_NAME,
                               SUBFOLDER_NAME=NULL,model=NULL){
  # --- Function to extract the input or selected predictors ---
  # Inputs: 
  #       - project_wd: string, path to the project folder
  #       - FOLDER_NAME: string, folder output from Cephalopod 
  #       - SUBFOLDER_NAME: string (optional), worms_id/scientificname of selected species
  #       - model: vector (optional), specific model(s) to retrieve predictors for
  # Outputs: 
  #       - predictors: tibble
  
  
  # 1. Retrieve the input predictors if there isn't a subfolder indicated
  if (is.null(SUBFOLDER_NAME)){
    load(paste0(project_wd, FOLDER_NAME,"/CALL.RData")) # load CALL
    predictors <- tibble(variable=CALL$ENV_VAR) # retrieve input environmental variables
  } 
  
  # 2. Retrieve the selected predictors for each model if there is a subfolder indicated
  else {
    load(paste0(project_wd, FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData")) 
    
    # models that were ran
    available_models <- names(MODEL)[!(names(MODEL) %in% c("MODEL_LIST","recommandations"))]
    
    # if model is NULL, use all available
    models_to_use <- if (is.null(model)) {
      available_models
    } else {
      intersect(available_models, model)
    }
    
    # get the predictor importance values for all the folds and bootstraps 
    predictors <- map_dfr(models_to_use, function(mod_name) {
      vip <- MODEL[[mod_name]]$vip %>% 
        # get the mean, median and std 
        group_by(variable) %>%
        summarise(
          mean = mean(value, na.rm = TRUE),
          median = median(value, na.rm = TRUE),
          sd = sd(value, na.rm = TRUE)) %>%
        # add the model name 
        mutate(model = mod_name)
    })
  }
  
  # 3. Return output 
  return(predictors)
}


# ----- extract_prediction -----

get_long_df <- function(arr, var_name, month_ids, factor_ids) {
  # --- Flatten an array ---
  # Input: 
  #     - arr: array to flatten 
  #     - var_name: name of the values in the array 
  # Output: 
  #     - array in long format 
  # grid names are necessary, just the row numbers for one of the month and varfactor works
  dimnames(arr) <- list(grid=seq(1,length(arr[,1,1])), month = month_ids, varfactor = factor_ids)
  as_tibble(as.table(arr)) %>%
    rename(!!var_name := n) 
}



extract_prediction <- function(project_wd,FOLDER_NAME,SUBFOLDER_NAME,
                               ensemble = TRUE, model=NULL,m=1:12,f=NULL){
  # --- Function to extract the projection, standard deviation and land mask ---
  # requires get_long_df function to be loaded 
  # Inputs: 
  #       - project_wd: string, path to the project folder
  #       - FOLDER_NAME: string, folder output from Cephalopod 
  #       - SUBFOLDER_NAME: string, worms_id/scientificname of selected species
  #       - ensemble: boolean, TRUE returns an ensemble model from the models that pass the QC (default) or indicated in model, FALSE returns a list of dataframes with all the models projections
  #       - model: vector (optional), specific model(s) to retrieve prediction for
  #       - m: vector (optional), months to pool average over e.g.: c(1,2,3), default to all
  #       - f: int (optional), one factor to plot (as a numeric) e.g.: 1 
  # Outputs: 
  #       - y_ens: data frame, with coordinates, and predicted mean, std, normalized standard deviation and coefficient of variation 
  
  # 1. Load the model output
  load(paste0(project_wd, FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # 2 Structure 
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # unpack the raster
  r0 <- CALL$ENV_DATA[[1]][[1]] # get the first environmental file to serve as structure
  r0_name <- CALL[["ENV_VAR"]][1]
  
  
  # 3. Model projections
  # models that were ran
  available_models <- names(MODEL)[!(names(MODEL) %in% c("MODEL_LIST","recommandations","ENSEMBLE"))]
  
  # if model is NULL, use all available (ensemble = FALSE) or the ones that pass the QC (ensemble = TRUE)
  models_to_use <- if (is.null(model)) {
    if (ensemble == FALSE){
      available_models
    } else {
      if (length(MODEL$MODEL_LIST) == 0){
        stop("No model passed the QC so we can't create a default ensemble. If you want a specific ensemble, indicate a vector of models with the model argument.")
      } else {
        MODEL$MODEL_LIST
      }
    }
  } else {
    intersect(available_models, model)
  }
  cat("Retrieving data for models:",models_to_use)
  
  
  # 3. Retrieve the data (beginning identical whatever ensemble)
  
  # loop over models 
  y_ens <- map_dfr(models_to_use, function(mod_name) {
    
    # message
    cat("\n Processing model:",mod_name)
    
    # extract array: [location, bootstrap, month, (varfactor)]
    y_hat <- MODEL[[mod_name]][["proj"]][["y_hat"]]
    
    # add 4th dimension if 3D array (no varfactor) => [location, bootstrap, month, 1]
    if (length(dim(y_hat)) == 3) {
      y_hat <- array(y_hat, dim = c(dim(y_hat), 1))
    }
    
    # determine months and factors to keep
    month_ids <- if (is.null(m)) seq_len(dim(y_hat)[3]) else m # all the dimensions or provided m
    factor_ids <- if (is.null(f)) seq_len(dim(y_hat)[4]) else f # all dimansions or provided f
    
    # average over bootstraps (dimension 2) => result is [location, month, varfactor]
    y_avg <- apply(y_hat[drop=FALSE], c(1, 3, 4), function(x)(x = mean(x, na.rm = TRUE)))  
    y_sd <- apply(y_hat[drop=FALSE], c(1, 3, 4), function(x)(x = sd(x, na.rm = TRUE)))  
    
    
    # 3.1. All models separately if ensemble = FALSE
    if (ensemble == FALSE) {
      
      # NSD: normalize sd by max of mean => [location,month,varfactor]
      model_max <- max(y_avg,na.rm=TRUE) 
      nsd <- y_sd/model_max
      cv <- y_sd/y_avg 
      cat("\n Normalized standard deviation for",mod_name, "was calculated with: ",model_max)
      
      # flatten all arrays: [grid, month, varfactor] → long format
      df <- reduce(
        list(
          get_long_df(y_avg[,month_ids,factor_ids,drop=FALSE], "predictedvalue", month_ids, factor_ids),
          get_long_df(y_sd[,month_ids,factor_ids,drop=FALSE],  "predictedsd", month_ids, factor_ids),
          get_long_df(nsd[,month_ids,factor_ids,drop=FALSE],   "predictednsd", month_ids, factor_ids),
          get_long_df(cv[,month_ids,factor_ids,drop=FALSE],    "predictedcv", month_ids, factor_ids)),
        full_join,
        by = c("grid", "month", "varfactor"))
      
      # add coordinates
      coords <- as.data.frame(xyFromCell(r0, 1:dim(y_avg)[1])) %>% 
        mutate(grid=row_number()) %>% 
        rename(longitude=x,latitude=y)
      df <- df %>%
        mutate(grid=as.integer(grid)) %>%
        suppressMessages(left_join(coords)) %>% 
        # remove the grid row numbers
        dplyr::select(-grid) %>%
        # add the model name and the model nsd max
        mutate(model = mod_name,
               max_nsd=model_max) %>% 
        # remove land values 
        drop_na(predictedvalue)
      
      return(df)
      
    } else {
      
      # 3.2. Average model if ensemble = TRUE 
      
      # 3.2.1. Average per model 
      
      cat("\n Processing ensemble model")
      
      # average over months and varfactors (dimension 2,3) => result is [location]
      mean_avg <- apply(y_avg[,month_ids,factor_ids], 1, mean, na.rm = TRUE)  
      mean_sd <- apply(y_sd[,month_ids,factor_ids], 1, mean, na.rm = TRUE)  
      
      # convert to dataframe
      coords <- as.data.frame(xyFromCell(r0, 1:length(mean_avg)))
      names(coords) <- c("longitude", "latitude")
      
      df <- tibble(
        latitude = coords$latitude,
        longitude = coords$longitude,
        predictedvalue  = mean_avg,
        predictedsd     = mean_sd,
        varfactor       = paste(factor_ids, collapse = ","),
        month           = paste(month_ids, collapse = ","))
      return(df)
    }
  })
  
  # 3.2.2. Average over models
  
  if (ensemble == TRUE) {
    
    # average over models
    y_ens <- y_ens %>% 
      group_by(latitude,longitude,varfactor,month) %>% 
      suppressMessages(summarize(across(everything(), ~mean(.x, na.rm = TRUE)))) %>% 
      # remove land values 
      drop_na(predictedvalue)
    
    # calculate the normalized standard deviation between models and the coefficient of variation
    global_max <- max(y_ens$predictedvalue,na.rm=TRUE)
    cat("\n Normalized standard deviation for the ensemble model was calculated with: ",global_max)
    y_ens <- y_ens %>% 
      mutate(predictednsd = predictedsd / global_max,
             predictedcv = ifelse(predictedvalue == 0, NA, predictedsd / predictedvalue)) %>% 
      mutate(model = "ensemble",
             max_nsd=global_max)
  }
  return(y_ens)
}

# ----- extract_land -----

extract_land <- function(project_wd,FOLDER_NAME){
  # --- Function to extract the land mask ---
  # Inputs: 
  #       - project_wd: string, path to the project folder
  #       - FOLDER_NAME: string, folder output from Cephalopod 
  # Outputs: 
  #       - land_df: data frame, with coordinates and NA for land values 
  
  # Load the environmental data 
  load(paste0(project_wd,FOLDER_NAME,"/CALL.RData"))
  
  # Retrieve the structure of the first environmental variable  
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x)) # unpack the raster
  r0 <- CALL$ENV_DATA[[1]][[1]] # get the first environmental file to serve as structure
  r0_name <- CALL[["ENV_VAR"]][1] # retrieve the name to update it 
  
  # Extract the land from the environmental file to use as landmask: 9999 for land and NA for ocean
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # Turn into a dateframe 
  land_df <- as.data.frame(land, xy = TRUE, na.rm = TRUE) %>% 
    rename(land=all_of(r0_name),
           longitude=x,
           latitude=y)
  
  return(land_df)
}
