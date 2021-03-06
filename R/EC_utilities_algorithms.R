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

    biomod_eval_method <- c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY",
                            "BIAS", "POD", "CSI", "ETS"), #vector of evaluation metrics

    # available from dismo::evaluate.R. Not originally implemented in biomod2::Evaluate.models.R
    dismo_eval_method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR"),

    rescale_all_models = a$rescale_all_models,  # if TRUE- model prediction scaled with binomial
    do_full_models = a$do_full_models,  # if TRUE, evaluation of whole dataset
    model_id = a$modeling_id,   # default= random; name of modeling procedure

    # matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which
    # part of data must be used for models calibration (TRUE) and for models
    # validation (FALSE). Each column correspund to a "RUN"
    # biomod.DataSplitTable = NULL

    species_name = response_info$occur_species,
    projection_name = "current",
    species_algo_str = ifelse(is.null(a$subset),
                              sprintf(model_algorithm,"%s_",
                                      response_info$occur_species),
                              sprintf(model_algorithm,"%_%s",
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


EC_compute_algorithm <- function( predictor_info,  # from EC_build_predictor
                                 pa_number_point,
                                 a,  # formerly EC.params
                                 EC_options_algorithm){  # from EC_options_algorithm
  # Define model options and compute the model
  # 1. Format the data as required by the biomod package
  model_data <- EC_format_biomod2 ( true.absen             = predictor_info$absen,
                                    pseudo.absen.points    = pa_number_point,
                                    pseudo.absen.strategy  = a$pa_strategy,
                                    pseudo.absen.disk.min  = a$pa_disk_min,
                                    pseudo.absen.disk.max  = a$pa_disk_max,
                                    pseudo.absen.sre.quant = a$pa_sre_quant,
                                    climate.data           = predictor_info$current_climate,
                                    occur                  = predictor_info$occur,
                                    species.name           = EC_options_algorithm$species_name,
                                    species_algo_str       = EC_options_algorithm$species_algo_str)

  # 2. Compute the model
  model_sdm <-
    biomod2::BIOMOD_Modeling (data               = model_data,
                              models             = model_algorithm,
                              models.options     = model_options,
                              NbRunEval          = EC_options_algorithm$NbRunEval,
                              DataSplit          = EC_options_algorithm$DataSplit,
                              Yweights           = EC_options_algorithm$Yweights,
                              Prevalence         = EC_options_algorithm$Prevalence,
                              VarImport          = EC_options_algorithm$VarImport,
                              models.eval.meth   = EC_options_algorithm$biomod_eval_meth,
                              SaveObj            = TRUE,
                              rescal.all.models  = EC_options_algorithm$rescal_all_models,
                              do.full.models     = EC_options_algorithm$do_full_models,
                              modeling.id        = EC_options_algorithm$model_id)
}


#_____________________________________________________________________________


save_model_object <- function ( model_compute,  # from model_compute
                                model_algorithm,
                                EC_options_algorithm){  # from build_biomod_functions
  # Save a variable importance plot (VIP) and the model object

  # save a VIP plot
  x.data <- attr(model_compute$model_data, "data.env.var")
  y.data <- attr(model_compute$model_data, "data.species")
  data1 <- data.frame(y.data, x.data)
  EC_plot_VIP (method = model_algorithm, data1 = data1, pdf = TRUE,
              filename = paste('vip_plot', EC_options_algorithm$species_algo_str,
                               sep = "_"),
              this.dir = paste(EC_options_algorithm$species_name,
                              "/models/EcoCommons", sep = ""))
  # save out the model object
  EC_save (model_compute$model_sdm, name = "model.object.RData")
}


#_____________________________________________________________________________


project_model_current_constrained <- function (predictor_info,  # from build_predictor
                                               constraint_info,
                                               model_compute,  # from model_compute
                                               EC_options_algorithm,  # from model_options_algorithm
                                               response_info,  # from build_response
                                               a,   # formerly EC.params
                                               model_algorithm){
  # Project over climate scenario. Constraint to continuous data layers

  if (predictor_info$genUnconstraintMap &&
      all(predictor_info$type == 'continuous') &&
      (!is.null(constraint_info$constraints) || predictor_info$generateCHull)) {
    model.proj <-
      biomod2::BIOMOD_Projection (modeling.output     = model_compute$model_sdm,
                                  new.env             = predictor_info$current_climate_orig,
                                  proj.name           = EC_options_algorithm$projection_name,
                                  xy.new.env          = EC_options_algorithm$xy_new_env,
                                  selected.models     = EC_options_algorithm$selected_models,
                                  binary.meth         = EC_options_algorithm$binary_method,
                                  filtered.meth       = EC_options_algorithm$filtered.method,
                                  #compress            = EC_options_algorithm$compress,
                                  build.clamping.mask = EC_options_algorithm$build_clamping_mask,
                                  silent              = EC_options_algorithm$silent,
                                  do.stack            = EC_options_algorithm$stack,
                                  keep.in.memory      = EC_options_algorithm$keep_in_memory,
                                  output.format       = EC_options_algorithm$output_format,
                                  on_0_1000           = FALSE)

    # remove the current climate rasters to release disk space
    EC_raster_remove (predictor_info$current_climate_orig)

    # convert projection output from grd to gtiff
    EC_GRD_to_GTIFF (file.path(getwd(),
                             response_info$occur_species,
                             paste("proj", EC_options_algorithm$projection_name,
                                   sep="_")),
                  
                   algorithm = ifelse(is.null(a$subset), model_algorithm,
                                      sprintf(model_algorithm, "_%s", a$subset)),
                  
                   filename_ext = "unconstrained")

    # save the projection
    EC_save_projection (model.proj, EC_options_algorithm$species_algo_str,
                        filename_ext="unconstrained")
  }
}


#_____________________________________________________________________________


project_model_current <- function( model_compute,  # from model_compute
                                   predictor_info,  # from build_predictor
                                   EC_options_algorithm,  # from model_options_algorithm
                                   response_info,  # from build_response
                                   a,
                                   model_algorithm ){
  # Project over climate scenario, without constraint

  model.proj <-
    biomod2::BIOMOD_Projection(modeling.output     = model_compute$model_sdm,
                               new.env             = predictor_info$current_climate,
                               proj.name           = EC_options_algorithm$projection_name,
                               xy.new.env          = EC_options_algorithm$xy_new_env,
                               selected.models     = EC_options_algorithm$selected_models,
                               binary.meth         = EC_options_algorithm$binary_method,
                               filtered.meth       = EC_options_algorithm$filtered.method,
                               #compress           = EC_options_algorithm$compress,
                               build.clamping.mask = EC_options_algorithm$build_clamping_mask,
                               silent              = EC_options_algorithm$silent,
                               do.stack            = EC_options_algorithm$stack,
                               keep.in.memory      = EC_options_algorithm$keep_in_memory,
                               output.format       = EC_options_algorithm$output_format,
                               on_0_1000           = FALSE)

  # remove the current_climate to release disk space
  EC_raster_remove (predictor_info$current_climate)

  # convert projection output from grd to gtiff
  EC_GRD_to_GTIFF (file.path(getwd(),
                           response_info$occur_species,
                           paste("proj", EC_options_algorithm$projection_name,
                                 sep="_")),
                 algorithm = ifelse(is.null(a$subset), model_algorithm,
                                    sprintf(model_algorithm,"_%s", a$subset)))


  # output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
  loaded.model <- biomod2::BIOMOD_LoadModels (model_compute$model_sdm,
                                              models = model_algorithm)

  EC_save_model_eval (loaded.model, model_compute$model_sdm,
                      EC_options_algorithm$species_algo_str)

  # save the projection
  EC_save_projection (model.proj, EC_options_algorithm$species_algo_str)

}


#_____________________________________________________________________________
## Subfunction within project_model_current_constrained() and project_model_current


EC_GRD_to_GTIFF <- function(folder,
                            algorithm,
                            filename_ext=NULL,
                            noDataValue=NULL) {
  # Convert all .gri/.grd found in folder to gtiff
  grdfiles <- list.files(path=folder,
                         pattern="^.*\\.grd")
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
    
    basename = paste(grdname, algorithm, sep="_")
    if (!is.null(filename_ext)) {
      basename = paste(grdname, algorithm, filename_ext, sep="_")
    }
    filename = file.path(folder, paste(basename, 'tif', sep="."))
    
    dtype = raster::dataType(grd)
    if (is.null(noDataValue)) {
      raster::writeRaster(grd, filename, datatype=dataType(grd),
                          format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
    }
    else {
      raster::writeRaster(grd, filename, datatype=dataType(grd), NAflag=noDataValue,
                          format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
    }
    file.remove(file.path(folder, paste(grdname, c("grd","gri"), sep=".")))  # remove grd files
  }
}