#' model_glm
#' 
#' Function to run a glm within EcoCommons
#' 
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly bccvl.params & bccvl.env)
# @param bccvl.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param bccvl.env Environment object; as for bccvl.params. Consider converting to S3
# @param occur.data filename for occurence data - generated from bccvl.params
# @param absen.data filename for absence data - generated from bccvl.params
# @param enviro.data.current list of filenames for climate data - generated from bccvl.params
# @param enviro.data.type continuous - generated from bccvl.params
# @param opt.tails predict parameter - not used anywhere in script
# @param outputdir root folder for output data - generated from bccvl.env
#'
#' @export model_glm


model_glm <- function(x){ # bccvl.params, bccvl.env){


  #define the working directory
  # NOTE: This section commented out in source code
  # directories <- set_directories(x$env)

  # set data for modelling
  response_info <- build_response_glm(x$params)
  predictor_info <- build_predictor_glm(x$params)

  # parameters to perform any biomod modelling
  model_options_glm <- build_glm_options(x$params)
  model_options_biomod <- build_biomod_options(x$params, response_info)

  ## MW CLEANED UP TO HERE ##

  # model accuracy statistics
  # these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
  dismo.eval.method = c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")

  # model accuracy statistics - combine stats from dismo and biomod2 for consistent output
  model.accuracy = c(dismo.eval.method, biomod.models.eval.meth)
  # TODO: these functions are used to evaluate the model ... configurable?

  # read current climate data
  current.climate.scenario = bccvl.enviro.stack(enviro.data.current, enviro.data.type, enviro.data.layer, resamplingflag=enviro.data.resampling)

  ###read in the necessary observation, background and environmental data
  occur = bccvl.species.read(occur.data, month.filter) #read in the observation data lon/lat
  absen = bccvl.species.read(absen.data, month.filter) #read in the observation data lon/lat

  # geographically constrained modelling
  if (!is.null(enviro.data.constraints) || enviro.data.generateCHall) {
    constrainedResults = bccvl.sdm.geoconstrained(current.climate.scenario, occur, absen, enviro.data.constraints, enviro.data.generateCHall);

    # Save a copy of the climate dataset
    current.climate.scenario.orig <- current.climate.scenario  
    current.climate.scenario <- constrainedResults$raster
    occur <- constrainedResults$occur
    absen <- constrainedResults$absen
  }

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio = bccvl.params$pa_ratio
  pa_number_point = 0
  if (pa_ratio > 0) {
    pa_number_point = floor(pa_ratio * nrow(occur))
  }

  ###run the models and store models
  ############### BIOMOD2 Models ###############
  # 1. Format the data
  # 2. Define the model options
  # 3. Compute the model
  # NOTE: Model evaluation is included as part of model creation

  # BIOMOD_Modeling(data, models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT'), models.options = NULL,
  #	NbRunEval=1, DataSplit=100, Yweights=NULL, Prevalence=NULL, VarImport=0, models.eval.meth = c('KAPPA','TSS','ROC'),
  #	SaveObj = TRUE, rescal.all.models = TRUE, do.full.models = TRUE, modeling.id = as.character(format(Sys.time(), '%s')),
  #	...)
  #
  # data	BIOMOD.formated.data object returned by BIOMOD_FormatingData
  # models vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF' and 'MAXENT'
  # models.options BIOMOD.models.options object returned by BIOMOD_ModelingOptions
  # NbRunEval	Number of Evaluation run
  # DataSplit	% of data used to calibrate the models, the remaining part will be used for testing
  # Yweights response points weights
  # Prevalence either NULL (default) or a 0-1 numeric used to build 'weighted response weights'
  # VarImport	Number of permutation to estimate variable importance
  # models.eval.meth vector of names of evaluation metric among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI' and 'ETS'
  # SaveObj keep all results and outputs on hard drive or not (NOTE: strongly recommended)
  # rescal.all.models	if true, all model prediction will be scaled with a binomial GLM
  # do.full.models if true, models calibrated and evaluated with the whole dataset are done
  # modeling.id character, the ID (=name) of modeling procedure. A random number by default.
  # ... further arguments :
  # DataSplitTable : a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a 'RUN'. If filled, args NbRunEval, DataSplit and do.full.models will be ignored.

  ###############
  #
  # GLM - generalized linear model (glm)
  #
  ###############

  # myBiomodOptions <- BIOMOD_ModelingOptions(GLM = list(type = 'quadratic', interaction.level = 0, myFormula = NULL,
  #	test = 'BIC', family = 'binomial', control = glm.control(epsilon = 1e-08, maxit = 1000, trace = FALSE)))
  # myFormula	: a typical formula object (see example). If not NULL, type and interaction.level args are switched off
  #	You can choose to either:
  #	1) generate automatically the GLM formula by using the type and interaction.level arguments
  #		type : formula given to the model ('simple', 'quadratic' or 'polynomial')
  #		interaction.level : integer corresponding to the interaction level between variables considered. Consider that
  #			interactions quickly enlarge the number of effective variables used into the GLM
  #	2) or construct specific formula
  # test : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.
  # family : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See family for details of family functions.) BIOMOD only runs on presence-absence data so far, so binomial family by default.
  # control : a list of parameters for controlling the fitting process. For glm.fit this is passed to glm.control
  #	glm.control(epsilon = 1e-8, maxit = 25, trace = FALSE)
  #		epsilon	- positive convergence tolerance e; the iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < e
  #		maxit - integer giving the maximal number of IWLS iterations
  #		trace - logical indicating if output should be produced for each iteration

  # 1. Format the data as required by the biomod package
  model.data = bccvl.biomod2.formatData(true.absen         = absen,
                                    pseudo.absen.points    = pa_number_point,
                                    pseudo.absen.strategy  = bccvl.params$pa_strategy,
                                    pseudo.absen.disk.min  = bccvl.params$pa_disk_min,
                                    pseudo.absen.disk.max  = bccvl.params$pa_disk_max,
                                    pseudo.absen.sre.quant = bccvl.params$pa_sre_quant,
                                    climate.data           = current.climate.scenario,
                                    occur                  = occur,
                                    species.name           = biomod.species.name,
                                    species_algo_str       = species_algo_str)

  # 2. Define the model options
  model.options <- BIOMOD_ModelingOptions(GLM = model.options.glm)
  # 3. Compute the model
  model.sdm <-
      BIOMOD_Modeling(data = model.data,
                      models=c('GLM'),
                      models.options=model.options,
                      NbRunEval=biomod.NbRunEval,
                      DataSplit=biomod.DataSplit,
                      Yweights=biomod.Yweights,
                      Prevalence=biomod.Prevalence,
                      VarImport=biomod.VarImport,
                      models.eval.meth=biomod.models.eval.meth,
                      SaveObj=TRUE,
                      rescal.all.models = biomod.rescal.all.models,
                      do.full.models = biomod.do.full.models,
                      modeling.id = biomod.modeling.id
                      )

  # save the VIP plot
  x.data <- attr(model.data,"data.env.var")
  y.data <- attr(model.data,"data.species")
  data1 = data.frame(y.data,x.data)
  bccvl.VIPplot(method="glm", data1=data1, pdf=TRUE, 
                filename=paste('vip_plot', species_algo_str, sep="_"), 
                this.dir=paste(biomod.species.name, "/models/bccvl", sep=""))

  # model output saved as part of BIOMOD_Modeling() # EMG not sure how to retrieve
  #save out the model object
  bccvl.save(model.sdm, name="model.object.RData")

  # Do projection over current climate scenario without constraint only if all env data layers are continuous.
  if (enviro.data.genUnconstraintMap &&
     all(enviro.data.type == 'continuous') && 
     (!is.null(enviro.data.constraints) || enviro.data.generateCHall)) {
      model.proj <-
          BIOMOD_Projection(modeling.output     = model.sdm,
                            new.env             = current.climate.scenario.orig,
                            proj.name           = projection.name,
                            xy.new.env          = biomod.xy.new.env,
                            selected.models     = biomod.selected.models,
                            binary.meth         = biomod.binary.meth,
                            filtered.meth       = biomod.filtered.meth,
                            #compress            = biomod.compress,
                            build.clamping.mask = biomod.build.clamping.mask,
                            silent              = opt.biomod.silent,
                            do.stack            = opt.biomod.do.stack,
                            keep.in.memory      = opt.biomod.keep.in.memory,
                            output.format       = opt.biomod.output.format,
                            on_0_1000           = FALSE)
      
      # remove the current climate rasters to release disk space
      bccvl.remove.rasterObject(current.climate.scenario.orig)

      # convert projection output from grd to gtiff
      bccvl.grdtogtiff(file.path(getwd(),
                                 biomod.species.name,
                                 paste("proj", projection.name, sep="_")), 
                       algorithm=ifelse(is.null(bccvl.params$subset), "glm", sprintf("glm_%s", bccvl.params$subset)),
                       filename_ext="unconstrained")

      # save the projection
      bccvl.saveProjection(model.proj, species_algo_str, filename_ext="unconstrained")
  }

  # predict for current climate scenario
  model.proj <-
      BIOMOD_Projection(modeling.output=model.sdm,
                        new.env=current.climate.scenario,
                        proj.name  = projection.name,  #basename(enviro.data.current), {{ species }}
                        xy.new.env = biomod.xy.new.env,
                        selected.models = biomod.selected.models,
                        binary.meth = biomod.binary.meth,
                        filtered.meth = biomod.filtered.meth,
                        #compress = biomod.compress,
                        build.clamping.mask = biomod.build.clamping.mask,
                        silent = opt.biomod.silent,
                        do.stack = opt.biomod.do.stack,
                        keep.in.memory = opt.biomod.keep.in.memory,
                        output.format = opt.biomod.output.format,
                        on_0_1000 = FALSE)

  # remove the current.climate.scenario to release disk space
  bccvl.remove.rasterObject(current.climate.scenario)

  # convert projection output from grd to gtiff
  bccvl.grdtogtiff(file.path(getwd(),
                             biomod.species.name,
                             paste("proj", projection.name, sep="_")),
                   algorithm=ifelse(is.null(bccvl.params$subset), "glm", sprintf("glm_%s", bccvl.params$subset)))


  # output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
  loaded.model = BIOMOD_LoadModels(model.sdm, models="GLM")
  bccvl.saveBIOMODModelEvaluation(loaded.model, model.sdm, species_algo_str) 	# save output

  # save the projection
  bccvl.saveProjection(model.proj, species_algo_str)

  return() # not set yet

} # end model_glm()



# subfunctions, listed in the order in which they appear within model_glm()

# NOTE: set_directories() currently not used
set_directories <- function(a){
  list(
    script = normalizePath(a$scriptdir),
    input =  normalizePath(a$inputdir),
    output =  normalizePath(a$outputdir)
  )
}

build_response_glm <- function(a){
  list(
    occur_data = a$species_occurrence_dataset$filename,
    occur_species = a$species_occurrence_dataset$species,
    month_filter = a$species_filter,
    #define the the lon/lat of the background / psuedo absence points to use -- 2 column matrix of longitude and latitude
    absen_data = a$species_absence_dataset$filename
  )
}

build_predictor_glm <- function(a){
  list(
    # NOTE: all entries previously prefaced with 'enviro.data.'
    #define the current enviro data to use
    current = lapply(a$environmental_datasets, function(x) x$filename),
    #type in terms of continuous or categorical
    type = lapply(a$environmental_datasets, function(x) x$type),
    #layer names for the current environmental layers used
    layer = lapply(a$environmental_datasets, function(x) x$layer),
    #geographic constraints.
    constraints = readLines(a$modelling_region$filename),
    #Indicate to generate and apply convex-hull polygon of occurrence dataset to constraint
    generateCHall = ifelse(is.null(a$generate_convexhull), FALSE, as.logical(bccvl.params$generate_convexhull)),
    #Indicate whether to generate unconstraint map or not. True by default
    genUnconstraintMap = ifelse(is.null(a$unconstraint_map), TRUE, as.logical(a$unconstraint_map)),
    # resampling (up / down scaling) if scale_down is TRUE, return 'lowest'
    resampling = ifelse(is.null(a$scale_down) ||
                                    as.logical(a$scale_down),
                                    'highest', 'lowest')
  )
}

build_glm_options <- function(a){
  list(
  	type = a$type,	#"simple", "quadratic" or "polynomial"; switched off if myFormula is not NULL
  	interaction_level = a$interaction_level, #integer corresponding to the interaction level between variables considered; switched off if myFormula is not NULL
  	# myFormula = NULL, #specific formula; if not NULL, type and interaction.level are args are switched off
  	test = a$test, #"AIC", "BIC" or "none"
  	family = a$family, #"binomial", "gaussian", "gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"
  	mustart = a$mustart, #starting values for the vector of means
  	control = list(
  		epsilon = a$control_epsilon, #positive convergence tolerance e
  		maxit = a$control_maxit, #integer giving the maximal number of IWLS iterations
  		trace = a$control_trace #logical indicating if output should be produced for each iteration
  	)
  )
}


#modeling.output #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
#new.env #a set of explanatory variables onto which models will be projected; must match variable names used to build the models
#proj.name #a character defining the projection name (a new folder will be created with this name)
# pseudo absences
build_biomod_options <- function(
  a, # formerly bccvl.params,
  response_info # from build_response_glm()
){
  list(
    NbRunEval = a$nb_run_eval,  # default 10; n-fold cross-validation; ignored if DataSplitTable is filled
    DataSplit = a$data_split, # default 100; % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
    Yweights = NULL, #response points weights
    Prevalence = a$prevalence, #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
    VarImport = a$var_import, # default 0; number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
    #EMG this parameter needs to be specified in order to get VariableImportance metrics during model evaluation
    models_eval_method = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS"), #vector of evaluation metrics
    rescale_all_models = a$rescale_all_models, #if true, all model prediction will be scaled with a binomial GLM
    do_full_models = a$do_full_models, #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
    modeling.id = a$modeling_id,  #character, the ID (=name) of modeling procedure. A random number by default
    # biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
    # EMG Need to test whether a NULL values counts as an argument
    species_name = response_info$occur_species, # used for various path and file name generation
    projection.name = "current",  #basename(enviro.data.current)
    species_algo_str = ifelse(is.null(a$subset), 
                              sprintf("%s_glm", response_info$occur_species), 
                              sprintf("%s_glm_%s", response_info$occur_species, a$subset)),
    PA_nb_rep = 0,
    PA_nb_absences = 0,
    xy_new_env = NULL, #optional coordinates of new.env data. Ignored if new.env is a rasterStack
    selected_models = a$selected_models, #'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (eg, = grep('_RF', getModelsBuiltModels(myBiomodModelOut)))
    # EMG If running one model at a time, this parameter becomes irrevelant
    binary_method = NULL, #a vector of a subset of models evaluation method computed in model creation
    filtered_method = NULL, #a vector of a subset of models evaluation method computed in model creation
    compress = a$compress, # default 'gzip'; compression format of objects stored on your hard drive. May be one of `xz', `gzip' or NULL
    build_clamping_mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
    options = list(
      silent = FALSE, #logical, if TRUE, console outputs are turned off
      stack = TRUE, #logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack
      keep_in_memory = TRUE, #logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object
      output.format = NULL #'.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class, output is .RData defining projections saving format (on hard drive)
    )
  )
}
