####
##
##  INPUT:
##
##  occur.data ... filename for occurence data
##  absen.data  ... filename for absence/background data
##  enviro.data.current ... list of filenames for climate data
##  enviro.data.type    ... continuous
##
##  outputdir ... root folder for output data

# The lon/lat of the observation records -- minimum of 2 column matrix of longitude and latitude
occur.data = bccvl.params$species_occurrence_dataset$filename
occur.species = bccvl.params$species_occurrence_dataset$species
month.filter = bccvl.params$species_filter
# The the lon/lat of the background / absence points -- 2 column matrix of longitude and latitude
absen.data = bccvl.params$species_absence_dataset$filename
# Optional bias file from user
bias.data = bccvl.params$bias_dataset$filename

# The current enviro data -- list of layers' filenames
enviro.data.current = lapply(bccvl.params$environmental_datasets,
                             function(x) {
                                fname = x$filename
                                return(fname)
                             }
)
#type in terms of continuous or categorical
enviro.data.type = lapply(bccvl.params$environmental_datasets, function(x) x$type)
#layer names for the current environmental layers used
enviro.data.layer = lapply(bccvl.params$environmental_datasets, function(x) x$layer)
# optional geographic constraints
enviro.data.constraints = bccvl.params$modelling_region
if (is.null(enviro.data.constraints) || enviro.data.constraints == '') {
  enviro.data.constraints = NULL
} else {
  enviro.data.constraints = readLines(bccvl.params$modelling_region$filename)
}
#Indicate to generate and apply convex-hull polygon of occurrence dataset to constraint
enviro.data.generateCHall = ifelse(is.null(bccvl.params$generate_convexhull), FALSE, as.logical(bccvl.params$generate_convexhull))
#Indicate whether to generate unconstraint map or not. True by default
enviro.data.genUnconstraintMap = ifelse(is.null(bccvl.params$unconstraint_map), TRUE, as.logical(bccvl.params$unconstraint_map))
# resampling (up / down scaling) if scale_down is TRUE, return 'lowest'
enviro.data.resampling = ifelse(is.null(bccvl.params$scale_down) ||
                                as.logical(bccvl.params$scale_down),
                                'highest', 'lowest')

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

# The number of background points
pa_number_point = bccvl.params$maximumbackground

# Generate background data if not supplied
if (is.null(absen)) {
    absen = generate_background_data(current.climate.scenario@layers[[1]], pa_number_point)
}

# Generate biased background data if a bias file is supplied
bg_biased = absen
if (!is.null(bias.data)) {
  bg_biased = generate_biased_background_data(absen, bias.data, pa_number_point)
}

#===================================================================
## Create data exploration outputs
## Extract predictor variable values for occurrence and background data
occ_coords <- subset(occur[,1:2]) # need to have a dataframe with only lat and lon to extract predictor variable values
occ_env <- raster::extract(current.climate.scenario, occ_coords)
bg_env <- raster::extract(current.climate.scenario, bg_biased)

# Save tables
occ_env_table <- data.frame(cbind(occur, occ_env))
bg_env_table <- data.frame(cbind(bg_biased, bg_env))
write.csv(occ_env_table, file.path(bccvl.env$outputdir,  paste0('occurrence_environmental_', occur.species, '_maxent.csv')), row.names = TRUE)
write.csv(bg_env_table, file.path(bccvl.env$outputdir,  paste0('background_environmental_', occur.species, '_maxent.csv')), row.names = TRUE)

## Run functions
included_layers = enviro.data.layer[which(enviro.data.type =='continuous')]
exploratoryPlots.contin.occ(occ_env, included_layers, outputdir=bccvl.env$outputdir)
exploratoryPlots.contin.bg(bg_env, included_layers, outputdir=bccvl.env$outputdir)
included_layers = enviro.data.layer[which(enviro.data.type !='continuous')]
exploratoryPlots.categ.occ(occ_env, included_layers, outputdir=bccvl.env$outputdir)
exploratoryPlots.categ.bg(bg_env, included_layers, outputdir=bccvl.env$outputdir)

#===================================================================

## Maxent model

# The Maxent function in the 'dismo' package is set up as follows:
# maxent(x, p, a=NULL, nbg=10000, factors=NULL, removeDuplicates=TRUE, ...)

# It uses the following arguments:
# x: Predictors, this is the RasterStack with the predictor layers,
# p: Occurrence data,
# a: Background points, 
# nbg: Number of background points. These are sampled randomly from the cells that are not NA in the first predictor variable. 
# This argument is ignored if background points are specified with argument a
# factors: which predictor variables, if any, should be considered as categorical?
# args: additonal arguments

# Set all Maxent arguments
maxent_args <- c(
  #duplicate records
  paste0('removeduplicates=', bccvl.params$removeduplicates), #remove duplicate presence records. If environmental data are in grids, duplicates are records in the same grid cell, otherwise, duplicates are records with identical coordinates.
  
  #background records
  paste0('maximumbackground=', bccvl.params$maximumbackground), #if the number of background points/grid cells is larger than this number, then this number of cells is chosen randomly for background points.
  paste0('addsamplestobackground=', bccvl.params$addsamplestobackground), #add to the background any sample for which has a combination of environmental values that isn't already present in the background
  paste0('addallsamplestobackground=', bccvl.params$addallsamplestobackground), #add all samples to the background, even if they have combinations of environmental values that are already present in the background
  
  #missing data
  paste0('allowpartialdata=', bccvl.params$allowpartialdata), #during model training, allow use of samples that have nodata values for one or more environmental variables
  
  #variable importance
  paste0('jackknife=', bccvl.params$jackknife), #default=FALSE; measure importance of each environmental variable by training with each environmental variable first omitted, then used in isolation.
  
  #random seed
  paste0('randomseed=', bccvl.params$randomseed), #if selected, a different random seed will be used for each run, so a different random test/train partition will be made and a different random subset of the background will be used, if applicable.
  
  #prevalence
  paste0('defaultprevalence=', bccvl.params$defaultprevalence), #default prevalence of the species: probability of presence at ordinary occurrence points. See Elith et al. Diversity and Distributions, 2011 for details
  
  #train/test settings
  paste0('randomtestpoints=', bccvl.params$randomtestpoints), #percentage of presence localities to be randomly set aside as test points, used to compute AUC, omission, etc.
  paste0('replicates=', bccvl.params$replicates), #number of replicate runs to do when cross-validating, bootstrapping or doing sampling with replacement runs.
  paste0('replicatetype=', bccvl.params$replicatetype), #if replicates > 1, do multiple runs of this type: 
  paste0('maximumiterations=', bccvl.params$maximumiterations), #stop training after this many iterations of the optimization algorithm
  paste0('convergencethreshold=', bccvl.params$convergencethreshold), #stop training when the drop in log loss per iteration drops below this number 
  
  #feature selection
  paste0('autofeature=', bccvl.params$autofeature), #automatically select which feature classes to use, based on number of training samples
  paste0('linear=', bccvl.params$linear), #allow linear features to be used
  paste0('quadratic=', bccvl.params$quadratic), #allow quadratic features to be used
  paste0('product=', bccvl.params$product), #allow product features to be used
  paste0('threshold=', bccvl.params$threshold), #allow threshold features to be used
  paste0('hinge=', bccvl.params$hinge), #allow hinge features to be used
  
  #feature settings
  paste0('lq2lqptthreshold=', bccvl.params$lq2lqptthreshold), #number of samples at which product and threshold features start being used
  paste0('l2lqthreshold=', bccvl.params$l2lqthreshold), #number of samples at which quadratic features start being used
  paste0('hingethreshold=', bccvl.params$hingethreshold), #number of samples at which hinge features start being used
  
  #regularization settings
  paste0('betamultiplier=', bccvl.params$betamultiplier), #multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.
  paste0('beta_threshold=', bccvl.params$beta_threshold), #regularization parameter to be applied to all threshold features; negative value enables automatic setting
  paste0('beta_categorical=', bccvl.params$beta_categorical), #regularization parameter to be applied to all categorical features; negative value enables automatic setting
  paste0('beta_lqp=', bccvl.params$beta_lqp), #regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
  paste0('beta_hinge=', bccvl.params$beta_hinge), #regularization parameter to be applied to all hinge features; negative value enables automatic setting
  
  #outputs - These are not shown in the UI, so unable to be changed by user
  'responsecurves=TRUE', #NB. default=FALSE; create graphs showing how predicted relative probability of occurrence depends on the value of each environmental variable
  'responsecurvesexponent=FALSE', #instead of showing the logistic value for the y axis in response curves, show the exponent (a linear combination of features).
  'pictures=TRUE', #create a .png image for each output grid
  'outputformat=raw', #representation of probabilities used in writing output grids, see Help for details
  'writeclampgrid=TRUE', #write a grid that shows the spatial distribution of clamping. At each point, the value is the absolute difference between prediction values with and without clamping.
  'writemess=TRUE', #a multidimensional environmental similarity surface (MESS) shows where novel climate conditions exist in the projection layers. The analysis shows both the degree of novelness and the variable that is most out of range at each point.
  'writeplotdata=FALSE', #write output files containing the data used to make response curves, for import into external plotting software.
  'outputgrids=TRUE', #write output grids. Turning this off when doing replicate runs causes only the summary grids (average, std, deviation, etc) to be written, not those for the individual runs.
  'plots=TRUE', #write various plots for inclusion in .html output
  'logfile=maxent.log', #file name to be used for writing debugging information about a run in output directory
  #'applythresholdrule=Fixed cumulative value 1', #apply a threshold rule, generating a binary outputgrid in addition to the regular prediction grid. Use the full name of the threshold rule in Maxent's html output as the argument. For example 'applyThresholdRule=Fixed cumulative value 1'.
  'logscale=TRUE', #if selected, all pictures of models will use a logarithmic scale for color-coding
  'writebackgroundpredictions=FALSE', #write .csv file with predictions at background points
  'fadebyclamping=FALSE', #reduce prediction at each point in projections by the difference between clamped and non-clamped output at that point

  #projection settings NB. These are not shown in the UI, so unable to be changed by user
  
  'extrapolate=TRUE', #predict to regions of environmental space outside the limits encountered during training
  'doclamp=TRUE' #apply clamping when projecting
  
)

#--------------------------------------------------------------------------------------
# Run model for the raw output
writeLines('Running model for raw output ...')
model_raw <- maxent(current.climate.scenario, occ_coords, bg_biased, path=file.path(bccvl.env$outputdir, "raw"), 
    args=maxent_args)


# Get entropy value used to transform the raw output to logistc and cloglog
lam <- parse_lambdas(model_raw)

# Project the raw model to model-fitting grids
pred.raw <- predict(model_raw, current.climate.scenario, args=c(outputformat="raw"))

# Transform the raw output to logistic/cloglog output -> this is used to produce the habitat suitability maps
pred.cloglog2 <- rmaxent::to_cloglog(pred.raw, from="raw", lam$entropy)
pred.logistic2 <- rmaxent::to_logistic(pred.raw, from="raw", lam$entropy)

#------------------------------------------------------------------------------------------------------------
## Maxent model results
#------------------------------------------------------------------------------------------------------------
## Habitat suitability map for each output format
fname = paste0("proj_current_", occur.species, "_maxent_cloglog")
writeRaster(pred.cloglog2, filename=file.path(bccvl.env$outputdir, fname), format = "GTiff", overwrite=TRUE)
png(filename=file.path(bccvl.env$outputdir, paste0(fname, ".png")), width=800, height=800)
title = "Habitat Suitability projection (cloglog)"
plot(pred.cloglog2, main=title, xlab='longitude', ylab='latitude')
dev.off()

fname = paste0("proj_current_", occur.species, "_maxent_logistic")
writeRaster(pred.logistic2, filename=file.path(bccvl.env$outputdir, fname), format = "GTiff", overwrite=TRUE)
png(filename=file.path(bccvl.env$outputdir, paste0(fname, ".png")), width=800, height=800)
title = "Habitat Suitability projection (logistic)"
plot(pred.logistic2, main=title, xlab='longitude', ylab='latitude')
dev.off()


## Response curves
# Generate response curves for cloglog and logistic output
# This script converts the output maps and figures to log and clog log (raw >log> cloglog)
# These are the maxent models for the logistic and cloglog outputs
writeLines('Running model for cloglog output ...')
model_rcCloglog <- maxent(current.climate.scenario, occ_coords, bg_biased, path=file.path(bccvl.env$outputdir, "cloglog"), 
                          args=maxent_args, outputformat="cloglog")
writeLines('Running model for logistic output ...')
model_rcLogistic <- maxent(current.climate.scenario, occ_coords, bg_biased, path=file.path(bccvl.env$outputdir, "logistic"), 
                          args=maxent_args, outputformat="logistic")

# Copy the roc file and omission file from cloglog results
tfile = file.path(bccvl.env$outputdir, "raw", "plots", "species_roc.png")
if (file.exists(tfile)) {
  file.copy(tfile, file.path(bccvl.env$outputdir, paste0("roc_", occur.species, "_maxent.png")))
}
tfile = file.path(bccvl.env$outputdir, "raw", "plots", "species_omission.png")
if (file.exists(tfile)) {
  file.copy(tfile, file.path(bccvl.env$outputdir, paste0("omcom_", occur.species, "_maxent.png")))
}

# Plot response curves
generate_response_curve(model_rcCloglog, names(current.climate.scenario), occur.species, "cloglog")
generate_response_curve(model_rcLogistic, names(current.climate.scenario), occur.species, "logistic")

## Maxent results table
maxentResults <- read.csv(file.path(bccvl.env$outputdir, "raw", "maxentResults.csv"))
maxentResults <- t(maxentResults)
head(maxentResults)
write.csv(maxentResults, file.path(bccvl.env$outputdir, paste0("maxent_results_", occur.species, ".csv")))

## Multivariate Environmental Similarity Surface (MESS)
# NB. Default is full=true however this returns a raster brick with all variables. Changed to full=FALSE to get only MESS result.
mess <- suppressWarnings(mess(current.climate.scenario, occ_env, full=FALSE))
writeRaster(x=mess, file.path(bccvl.env$outputdir, paste0("mess_", occur.species, "_maxent.tif")),
            format="GTiff", overwrite=TRUE)
png(filename=file.path(bccvl.env$outputdir, paste0("mess_", occur.species, "_maxent.png")), 
    width=800, height=800)
plot(mess, col=terrain.colors(10, alpha = 1), xlab='longitude', ylab='latitude')
dev.off()

## Limiting Factors Map which is categorical data
limiting_factors <- limiting(current.climate.scenario, model_raw)
png(filename=file.path(bccvl.env$outputdir, paste0("limitfactors_", occur.species, "_maxent.png")),
    width=800, height=800)
lf_plot = levelplot(limiting_factors)
print(lf_plot)
dev.off()

#------------------------------------------------------------------------------------------------------------
## Maxent model evaluation
#------------------------------------------------------------------------------------------------------------
## Random NULL model

# Usage:
# x = data frame with environmental predictor values for occurrences
# model = model function that creates a model of class 'DistModel', example uses bioclim
# n = sample size
# rep = number of repetitions
# pa = boolean (TRUE or FALSE), presence-only or presence-background

writeLines('Random NULL model ...')
nr <- nullRandom(as.data.frame(occ_env), maxent, n=min(25, nrow(occ_env)-1), rep=bccvl.params$randomnullrep, pa=TRUE)

# Determine the 95% confidence interval (CI)
CI <- nr[[round(0.95*bccvl.params$randomnullrep)]]
print(CI)

# save the model AUC and null model AUC
auc_df = data.frame('Training Model AUC'=maxentResults['Training.AUC',], 'NULL MODEL AUC'=CI@auc)
write.csv(auc_df, file.path(bccvl.env$outputdir, paste0("nullmodel_", occur.species, "_maxent.csv")))

#-------------------------------------------------------------------------------
##  Evaluation statistics
equal_sens_spec_threshold = saveEvaluationStatistics(occ_coords, bg_biased, model_raw, current.climate.scenario, occur.species)
writeLines(paste0('equal_sens_spec_threshold = ', equal_sens_spec_threshold))

## Area/Extent of Occuppancy (AOO/EOO)
#computeAreaExtentOccuopancy(occur, occur.species)

# generate species threshold map based on cloglog projection 
proj_file = file.path(bccvl.env$outputdir, paste0("proj_current_", occur.species, "_maxent_cloglog.tif"))
proj_raster = bccvl.raster.load(proj_file)
bccvl.generateThresholdMap(proj_raster, equal_sens_spec_threshold)

#--------------------------------------------------------------------------------
# Do projection over current climate scenario without constraint only if all env data layers are continuous.
if (enviro.data.genUnconstraintMap && 
   all(enviro.data.type == 'continuous') && 
   (!is.null(enviro.data.constraints) || enviro.data.generateCHall)) {

   writeLines("Doing projection over unconstraint climate scenario ...")
   pred.raw.unconstraint <- predict(model_raw, current.climate.scenario.orig, args=c(outputformat="raw"))

   # Transform the raw output to logistic/cloglog output -> this is used to produce the habitat suitability maps
   pred.cloglog2_unconstraint <- rmaxent::to_cloglog(pred.raw.unconstraint, from="raw", lam$entropy)
   pred.logistic2_unconstraint <- rmaxent::to_logistic(pred.raw.unconstraint, from="raw", lam$entropy)

  ## Habitat suitability map for each output format
  fname = paste0("proj_current_", occur.species, "_maxent_cloglog_unconstraint")
  writeRaster(pred.cloglog2_unconstraint, filename=file.path(bccvl.env$outputdir, fname), format = "GTiff", overwrite=TRUE)
  fname = paste0("proj_current_", occur.species, "_maxent_logistic_unconstraint")
  writeRaster(pred.logistic2_unconstraint, filename=file.path(bccvl.env$outputdir, fname), format = "GTiff", overwrite=TRUE)
}

writeLines('Done')