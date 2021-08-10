#===================================================================
# subfunctions, for the Statistical Regression algorithms scripts

build_biomod_options <- function(
  a, # formerly EC.params
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



model_compute <- function(
  predictor_info,  # from build_predictor
  pa_number_point,
  a,  # formerly EC.params
  model_options_biomod,  # from build_biomod_functions
  model_options_glm ){  # from build_glm_options

  # Define model options and compute the model

  # 1. Format the data as required by the biomod package
  model_data <- EC_FormatDataBIOMOD2 (true.absen             = predictor_info$absen,
                                      pseudo.absen.points    = pa_number_point,
                                      pseudo.absen.strategy  = a$pa_strategy,
                                      pseudo.absen.disk.min  = a$pa_disk_min,
                                      pseudo.absen.disk.max  = a$pa_disk_max,
                                      pseudo.absen.sre.quant = a$pa_sre_quant,
                                      climate.data           = predictor_info$current_climate,
                                      occur                  = predictor_info$occur,
                                      species.name           = model_options_biomod$species_name,
                                      species_algo_str       = model_options_biomod$species_algo_str)

  # 2. Define the model options
  model_options <- BIOMOD_ModelingOptions(GLM = model_options_glm)

  # 3. Compute the model
  model_sdm <-
    BIOMOD_Modeling(data               = model_data,
                    models             = model_algorithm,
                    models.options     = model_options,
                    NbRunEval          = model_options_biomod$NbRunEval,
                    DataSplit          = model_options_biomod$DataSplit,
                    Yweights           = model_options_biomod$Yweights,
                    Prevalence         = model_options_biomod$Prevalence,
                    VarImport          = model_options_biomod$VarImport,
                    models.eval.meth   = model_options_biomod$biomod_eval_meth,
                    SaveObj            = TRUE,
                    rescal.all.models  = model_options_biomod$rescal_all_models,
                    do.full.models     = model_options_biomod$do_full_models,
                    modeling.id        = model_options_biomod$model_id)
}



save_model_object<-function(
  model_compute_sdm,  # from model_compute
  model_algorithm,
  model_options_biomod){  # from build_biomod_functions

  # Save a variable importance plot (VIP) and the model object

  # save a VIP plot
  x.data <- attr(model_compute_sdm$model_data, "data.env.var")
  y.data <- attr(model_compute_sdm$model_data, "data.species")
  data1 <- data.frame(y.data, x.data)
  EC_VIPplot(method = model_algorithm, data1 = data1, pdf = TRUE,
             filename = paste('vip_plot', model_options_biomod$species_algo_str,
                              sep = "_"),
             this.dir = paste(model_options_biomod$species_name,
                              "/models/EcoCommons", sep = ""))
  # save out the model object
  EC_Save(model_compute_sdm$model_sdm, name = "model.object.RData")
}



project_model_current_constrained <- function(
  predictor_info,  # from build_predictor
  model_compute_sdm,  # from model_compute
  model_options_biomod,  # from build_biomod_options
  response_info,  # from build_response
  a,   # formerly EC.params
  model_algorithm){

  # Project over climate scenario. Constraint to continuous data layers

  if (predictor_info$genUnconstraintMap &&
      all(enviro.data.type == 'continuous') &&
      (!is.null(constraints) || predictor_info$generateCHull)) {
    model.proj <-
      BIOMOD_Projection(modeling.output     = model_compute_sdm$model_sdm,
                        new.env             = predictor_info$current_climate_orig,
                        proj.name           = model_options_biomod$projection_name,
                        xy.new.env          = model_options_biomod$xy_new_env,
                        selected.models     = model_options_biomod$selected_models,
                        binary.meth         = model_options_biomod$binary_method,
                        filtered.meth       = model_options_biomod$filtered.method,
                        #compress            = model_options_biomod$compress,
                        build.clamping.mask = model_options_biomod$build_clamping_mask,
                        silent              = model_options_biomod$silent,
                        do.stack            = model_options_biomod$stack,
                        keep.in.memory      = model_options_biomod$keep_in_memory,
                        output.format       = model_options_biomod$output_format,
                        on_0_1000           = FALSE)

    # remove the current climate rasters to release disk space
    EC_RevRasterObject(predictor_info$current_climate_orig)

    # convert projection output from grd to gtiff
    EC_GRDtoGTIFF(file.path(getwd(),
                            response_info$occur_species,
                            paste("proj", model_options_biomod$projection_name,
                                  sep="_")),
                  algorithm = ifelse(is.null(a$subset), model_algorithm,
                                     sprintf(model_algorithm, "_%s", a$subset)),
                  filename_ext = "unconstrained")

    # save the projection
    EC_SaveProjection(model.proj, model_options_biomod$species_algo_str,
                      filename_ext="unconstrained")
  }
}



project_model_current<-function(
  model_compute_sdm,  # from model_compute
  predictor_info,  # from build_predictor
  model_options_biomod,  # from build_biomod_options
  response_info,  # from build_response
  a,
  model_algorithm ){

  # Project over climate scenario, without constraint

  model.proj <-
    BIOMOD_Projection(modeling.output     = model_compute_sdm$model_sdm,
                      new.env             = predictor_info$current_climate,
                      proj.name           = model_options_biomod$projection_name,
                      xy.new.env          = model_options_biomod$xy_new_env,
                      selected.models     = model_options_biomod$selected_models,
                      binary.meth         = model_options_biomod$binary_method,
                      filtered.meth       = model_options_biomod$filtered.method,
                      #compress           = model_options_biomod$compress,
                      build.clamping.mask = model_options_biomod$build_clamping_mask,
                      silent              = model_options_biomod$silent,
                      do.stack            = model_options_biomod$stack,
                      keep.in.memory      = model_options_biomod$keep_in_memory,
                      output.format       = model_options_biomod$output_format,
                      on_0_1000           = FALSE)

  # remove the current_climate to release disk space
  EC_RevRasterObject(predictor_info$current_climate)

  # convert projection output from grd to gtiff
  EC_GRDtoGTIFF(file.path(getwd(),
                          response_info$occur_species,
                          paste("proj", model_options_biomod$projection_name,
                                sep="_")),
                algorithm=ifelse(is.null(a$subset), model_algorithm,
                                 sprintf(model_algorithm,"_%s", a$subset)))


  # output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
  loaded.model <- BIOMOD_LoadModels (model_compute_sdm$model_sdm,
                                     models = model_algorithm)

  EC_SaveModelEval(loaded.model, model_compute_sdm$model_sdm,
                   model_options_biomod$species_algo_str)

  # save the projection
  EC_SaveProjection(model.proj, model_options_biomod$species_algo_str)

}
