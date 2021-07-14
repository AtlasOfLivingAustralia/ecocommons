############################################################
###     EcoCommons script to FDA - using biomod2     ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs the maximum entropy modelling (MaxEnt) algorithm 
##  (Philips et al. 2004) for Species Distribution Modelling, using the 
##  biomod2 package on R and the input data generated on the EcoCommons 
##  platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## BIOMOD2 Models
# general parameters to perform any biomod modelling

biomod.NbRunEval <- EC.params$nb_run_eval  # default 10; n-fold cross-validation; ignored if DataSplitTable is filled
biomod.DataSplit <- EC.params$data_split # default 100; % for calibrating/training, remainder for testing; ignored if DataSplitTable is filled
biomod.Yweights <- NULL #response points weights
biomod.Prevalence <- EC.params$prevalence #either NULL (default) or a 0-1 numeric used to build "weighted response weights"
biomod.VarImport <- EC.params$var_import # default 0; number of resampling of each explanatory variable to measure the relative importance of each variable for each selected model
#EMG this parameter needs to be specified in order to get VariableImportance metrics during model evaluation
biomod.models.eval.meth <- c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY",
                             "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics
biomod.rescal.all.models <- EC.params$rescale_all_models #if true, all model prediction will be scaled with a binomial GLM
biomod.do.full.models <- EC.params$do_full_models #if true, models calibrated and evaluated with the whole dataset are done; ignored if DataSplitTable is filled
biomod.modeling.id <- EC.params$modeling_id #character, the ID (=name) of modeling procedure. A random number by default
# biomod.DataSplitTable = NULL #a matrix, data.frame or a 3D array filled with TRUE/FALSE to specify which part of data must be used for models calibration (TRUE) and for models validation (FALSE). Each column correspund to a "RUN". If filled, args NbRunEval, DataSplit and do.full.models will be ignored
# EMG Need to test whether a NULL values counts as an argument
biomod.species.name <- occur.species # used for various path and file name generation
projection.name <- "current"  #basename(enviro.data.current)
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_fda", occur.species), 
                          sprintf("%s_fda_%s", occur.species,
                                  EC.params$subset))


# model-specific arguments to create a biomod model
model.options.fda <- list(
	method = EC.params$method #regression method used in optimal scaling; "polyreg", "mars", "bruto" or "gen.ridge"
)

############### BIOMOD2 Models ###############
#
# general parameters to project any biomod modelling
#
#modeling.output #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
#new.env #a set of explanatory variables onto which models will be projected; must match variable names used to build the models
#proj.name #a character defining the projection name (a new folder will be created with this name)
# pseudo absences
biomod.PA.nb.rep = 0
biomod.PA.nb.absences = 0

biomod.xy.new.env = NULL #optional coordinates of new.env data. Ignored if new.env is a rasterStack
biomod.selected.models = EC.params$selected_models #'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (eg, = grep('_RF', getModelsBuiltModels(myBiomodModelOut)))
# EMG If running one model at a time, this parameter becomes irrevelant
biomod.binary.meth = NULL #a vector of a subset of models evaluation method computed in model creation
biomod.filtered.meth = NULL #a vector of a subset of models evaluation method computed in model creation
biomod.compress = EC.params$compress # default 'gzip'; compression format of objects stored on your hard drive. May be one of `xz', `gzip' or NULL
biomod.build.clamping.mask = FALSE #if TRUE, a clamping mask will be saved on hard drive
opt.biomod.silent = FALSE #logical, if TRUE, console outputs are turned off
opt.biomod.do.stack = TRUE #logical, if TRUE, attempt to save all projections in a unique object i.e RasterStack
opt.biomod.keep.in.memory = TRUE #logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object
opt.biomod.output.format = NULL #'.Rdata', '.grd' or '.img'; if NULL, and new.env is not a Raster class, output is .RData defining projections saving format (on hard drive)


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy <- c(dismo.eval.method, biomod.models.eval.meth)


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
# FDA - flexible discriminant analysis (fda)
#
###############

# myBiomodOptions <- BIOMOD_ModelingOptions(FDA = list(method = 'mars'))
# method : regression method used in optimal scaling.
# Default is linear regression via the function polyreg, resulting in linear discriminant analysis.
# Other possibilities are mars and bruto. For Penalized Discriminant analysis gen.ridge is appropriate.

# 1. Format the data as required by the biomod package
model.data <- EC_FormatDataBIOMOD2(true.absen         = absen,
                                  pseudo.absen.points    = pa_number_point,
                                  pseudo.absen.strategy  = EC.params$pa_strategy,
                                  pseudo.absen.disk.min  = EC.params$pa_disk_min,
                                  pseudo.absen.disk.max  = EC.params$pa_disk_max,
                                  pseudo.absen.sre.quant = EC.params$pa_sre_quant,
                                  climate.data           = current.climate.scenario,
                                  occur                  = occur,
                                  species.name           = biomod.species.name,
                                  species_algo_str       = species_algo_str)

# 2. Define the model options
model.options <- BIOMOD_ModelingOptions(FDA = model.options.fda)
# 3. Compute the model
model.sdm <-
    BIOMOD_Modeling(data              = model.data,
                    models            = c('FDA'),
                    models.options    = model.options,
                    NbRunEval         = biomod.NbRunEval,
                    DataSplit         = biomod.DataSplit,
                    Yweights          = biomod.Yweights,
                    Prevalence        = biomod.Prevalence,
                    VarImport         = biomod.VarImport,
                    models.eval.meth  = biomod.models.eval.meth,
                    SaveObj           = TRUE,
                    rescal.all.models = biomod.rescal.all.models,
                    do.full.models    = biomod.do.full.models,
                    modeling.id       = biomod.modeling.id
                    )

#save out the model object
EC_Save(model.sdm, name="model.object.RData")

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

    # remove the current.climate.scenario to release disk space
    EC_RevRasterObject(current.climate.scenario.orig)
        
    # convert projection output from grd to gtiff
    EC_GRDtoGTIFF(file.path(getwd(),
                               biomod.species.name,
                               paste("proj", projection.name, sep="_")), 
                     algorithm=ifelse(is.null(EC.params$subset), "fda", sprintf("fda_%s", EC.params$subset)),
                     filename_ext="unconstrained")

    # save the projection
    EC_SaveProjection(model.proj, species_algo_str, filename_ext="unconstrained")
}

# predict for current climate scenario
model.proj <-
    BIOMOD_Projection(modeling.output     = model.sdm,
                      new.env             = current.climate.scenario,
                      proj.name           = projection.name,
                      xy.new.env          = biomod.xy.new.env,
                      selected.models     = biomod.selected.models,
                      binary.meth         = biomod.binary.meth,
                      filtered.meth       = biomod.filtered.meth,
                      #compress           = biomod.compress,
                      build.clamping.mask = biomod.build.clamping.mask,
                      silent              = opt.biomod.silent,
                      do.stack            = opt.biomod.do.stack,
                      keep.in.memory      = opt.biomod.keep.in.memory,
                      output.format       = opt.biomod.output.format,
                      on_0_1000           = FALSE)

# remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# convert projection output from grd to gtiff
EC_GRDtoGTIFF(file.path(getwd(),
                           biomod.species.name,
                           paste("proj", projection.name, sep="_")),
                 algorithm=ifelse(is.null(EC.params$subset), "fda", sprintf("fda_%s", EC.params$subset)))

# output is saved as part of the projection, format specified in arg 'opt.biomod.output.format'
# evaluate model
loaded.model <- BIOMOD_LoadModels(model.sdm, models="FDA") # load model
EC_SaveBIOMODModelEval(loaded.model, model.sdm, species_algo_str)

# save the projection
EC_SaveProjection(model.proj, species_algo_str)
