#' Subfunctions for all Species Distribution Modelling algorithms EXCEPT the 
#' Geographical ones. It includes build evaluation parameters, define and compute
#' the model following biomod2 format, saving and export functions
#' 
#' @importFrom biomod2 BIOMOD_LoadModels
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom raster dataType
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' 
#' @keywords internal

EC_options_algorithm <- function ( a, # formerly EC.params
                                  response_info, # from build_response()
                                  model_algorithm){
  
  # Build evaluation parameters for biomod2 and dismo modelling
  
  list(
    NbRunEval  = a$nb_run_eval,  # default 10; n-fold cross-validation
    DataSplit  = a$data_split,  # default 100; % for calibrating/training
    Yweights   = NULL,  # response points weights
    Prevalence = a$prevalence,  # default NULL; or 0-1 numeric used to build "weighted response weights"

    # number of resampling of each variable; measure the importance of each variable for each model
    VarImport  = a$var_import, # default 0

    # Model accuracy statistics
    # some are available from biomod2::Evaluate.models.R
    biomod_eval_method = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS"), #vector of evaluation metrics
    # others from dismo::evaluate.R 
    dismo_eval_method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR"),
    
    # model accuracy statistics - combine stats from dismo and biomod2 for consistent output
    model_accuracy = c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", "BIAS", 
                       "POD", "CSI", "ETS", "ODP", "TNR", "FPR", "FNR", "NPP",
                       "MCR", "OR"),
 
    rescale_all_models = a$rescale_all_models,  # if TRUE- model prediction scaled with binomial
    do_full_models = a$do_full_models,  # if TRUE, evaluation of whole dataset
    model_id = a$modeling_id,   # default= random; name of modeling procedure

    species_name = response_info$occur_species,
    projection_name = "current",
    species_algo_str = ifelse(is.null(a$subset),
                              paste0(model_algorithm,"_",
                                      response_info$occur_species),
                              paste0(model_algorithm,"_",
                                      response_info$occur_species, a$subset)),
    
    PA_nb_rep = 0,
    PA_nb_absences = 0,
    xy_new_env = NULL,  # optional coordinates. Ignored if new.env is a rasterStack
    selected_models = a$selected_models,  # default = all;  when a subset vector of modeling.output models  is to be computed
    binary_method = NULL,  # vector/subset of models evaluation method computed in model creation
    filtered_method = NULL,  # a vector of a subset of models evaluation method computed in model creation
    compress = a$compress,  # default 'gzip'; compress objects stored on your hard drive
    build_clamping_mask = FALSE,  # if TRUE, a clamping mask will be saved on hard drive
    options = list(
      silent = FALSE,  # logical; if TRUE, console outputs are turned off
      stack = TRUE,  # logical; if TRUE, attempt to save all projections in a unique object i.e RasterStack
      keep_in_memory = TRUE,  # logical; if FALSE only the link pointing is stored in output object
      output_format = NULL  # '.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class
    )
  )
}

#_____________________________________________________________________________


project_model_current_constrained <- function ( constraint_info, # from build_constraint
                                                predictor_info,  # from build_predictor
                                                model_sdm,       # from model created using biomod2::BIOMOD_Modeling
                                                dataset_info,    # from build_dataset
                                                model_options_algorithm, # from EC_options_algorithm
                                                model_algorithm,
                                                a) {             # formerly EC.params
                                                
  # Project over climate scenario. Constraint to continuous data layers

  if (constraint_info$genUnconstraintMap &&
      all(predictor_info$type == 'continuous') &&
      (!is.null(constraint_info$constraints) || constraint_info$generateCHull)) {
    model_proj <-
      biomod2::BIOMOD_Projection (modeling.output     = model_sdm,
                                  new.env             = dataset_info$current_climate_orig,
                                  proj.name           = model_options_algorithm$projection_name,
                                  xy.new.env          = model_options_algorithm$xy_new_env,
                                  selected.models     = model_options_algorithm$selected_models,
                                  binary.meth         = model_options_algorithm$binary_method,
                                  filtered.meth       = model_options_algorithm$filtered.method,
                                  build.clamping.mask = model_options_algorithm$build_clamping_mask,
                                  silent              = model_options_algorithm$silent,
                                  do.stack            = model_options_algorithm$stack,
                                  keep.in.memory      = model_options_algorithm$keep_in_memory,
                                  output.format       = model_options_algorithm$output_format,
                                  on_0_1000           = FALSE)

    # remove the current climate rasters to release disk space
    EC_raster_remove (dataset_info$current_climate_orig)

    # convert projection output from grd to gtiff
    EC_GRD_to_GTIFF (file.path(EC.env$outputdir,
                               dataset_info$occur_species,
                               paste("proj", model_options_algorithm$projection_name,
                                     sep="_")),
                  
                   algorithm = ifelse(is.null(a$subset), model_algorithm,
                                      sprintf(model_algorithm, "_%s", a$subset)),
                   filename_ext = "unconstrained")

    # save the projection
    EC_save_projection (model_proj,
                        model_options_algorithm$species_algo_str,
                        filename_ext = "unconstrained")
  }
}


#_____________________________________________________________________________


project_model_current <- function( model_sdm,               # from model created using biomod2::BIOMOD_Modeling
                                   dataset_info,            # from build_dataset
                                   model_options_algorithm, # from EC_options_algorithm
                                   model_algorithm) {       # formerly EC.params

  # Project over climate scenario, without constraint

  model_proj <-
    biomod2::BIOMOD_Projection(modeling.output     = model_sdm,
                               new.env             = dataset_info$current_climate,
                               proj.name           = model_options_algorithm$projection_name,
                               xy.new.env          = model_options_algorithm$xy_new_env,
                               selected.models     = model_options_algorithm$selected_models,
                               binary.meth         = model_options_algorithm$binary_method,
                               filtered.meth       = model_options_algorithm$filtered.method,
                               build.clamping.mask = model_options_algorithm$build_clamping_mask,
                               silent              = model_options_algorithm$silent,
                               do.stack            = model_options_algorithm$stack,
                               keep.in.memory      = model_options_algorithm$keep_in_memory,
                               output.format       = model_options_algorithm$output_format,
                               on_0_1000           = FALSE)

  # remove the current_climate to release disk space
  EC_raster_remove (dataset_info$current_climate)

  # convert projection output from grd to gtiff
  EC_GRD_to_GTIFF (file.path(EC.env$outputdir,
                           dataset_info$occur_species,
                           paste("proj", model_options_algorithm$projection_name,
                                 sep="_")),
                 algorithm = ifelse(is.null(a$subset),
                                    paste0(model_algorithm),
                                    paste0(model_algorithm,"_",
                                           a$subset)))
  
  # output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
  
  # extract the root of filenames used by biomod2 to save model results
  filenames <- dir(this.dir)
  
  
  loaded_model <- biomod2::BIOMOD_LoadModels (model_sdm,
                                              models = model_algorithm)

  EC_save_BIOMOD_model_eval (loaded_model,
                             model_sdm,
                             model_options_algorithm$species_algo_str)

  # save the projection
  EC_save_projection (model_proj,
                      model_options_algorithm$species_algo_str)

}


#_____________________________________________________________________________
## Subfunction within project_model_current_constrained() and project_model_current


EC_GRD_to_GTIFF <- function(folder,
                            algorithm,
                            filename_ext = NULL,
                            noDataValue  = NULL) {
  # Convert all .gri/.grd found in folder to gtiff
  
  grdfiles <- list.files(path    = folder,
                         pattern = "^.*\\.grd")
  for (grdfile in grdfiles) {
    
    ext <- file_ext(grdfile)
    if (!is.null(ext)) {
      pattern = paste0('\\.', ext, '$')
      grdname <- sub(pattern, '', grdfile)
    }
    
    grd <- raster::raster(file.path(folder, grdfile))
    
    if (is.na(proj4string(grd))) {
      crs = CRS("+init=epsg:4326")  # if projection is missing, initialise it to EPSG:4326
      proj4string(grd) <- crs
    }
    
    basename = paste(grdname, algorithm, sep = "_")
    if (!is.null(filename_ext)) {
      basename = paste(grdname, algorithm, filename_ext, sep="_")
    }
    filename = file.path(folder, paste(basename, 'tif', sep="."))
    
    dtype = raster::dataType(grd)
    
    if (is.null(noDataValue)) {
      raster::writeRaster(grd, filename, datatype = dataType(grd),
                          format = "GTiff", options = c("COMPRESS=LZW", "TILED=YES"),
                          overwrite = TRUE)
    }
    else {
      raster::writeRaster(grd, filename, datatype = dataType(grd), NAflag = noDataValue,
                          format = "GTiff", options = c("COMPRESS= LZW", "TILED= YES"),
                          overwrite = TRUE)
    }
    file.remove(file.path(folder, paste(grdname, c("grd","gri"), sep = ".")))  # remove grd files
  }
}