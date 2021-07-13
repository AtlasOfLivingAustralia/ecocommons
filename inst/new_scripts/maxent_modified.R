############################################################
###     EcoCommons script to MaxEnt - using biomod2     ###
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
## Load and edit your dataset with parameters specific to MaxEnt

# Set number of background points
pa_number_point <- EC.params$maximumbackground

# Generate background data, if not supplied
if (is.null(absen)) {
    absen <- EC_CreateBackgroundData(current.climate.scenario@layers[[1]], 
                                     pa_number_point)
}

# Generate biased background data if a bias file is supplied
bg_biased <- absen
if (!is.null(bias.data)) {
  bg_biased <- EC_CreateeBiasedBackgroundData(absen, bias.data, pa_number_point)
}

#===================================================================
## Create data exploration outputs

# Extract predictor variable values for occurrence and background data
occ_coords <- subset(occur[,1:2]) # data frame with only lat and lon to extract predictor values
occ_env <- raster::extract(current.climate.scenario, occ_coords)
bg_env <- raster::extract(current.climate.scenario, bg_biased)

# Save tables
occ_env_table <- data.frame(cbind(occur, occ_env))
bg_env_table <- data.frame(cbind(bg_biased, bg_env))
write.csv(occ_env_table, file.path(EC.env$outputdir,  
                                   paste0('occurrence_environmental_', 
                                        occur.species, '_maxent.csv')), row.names = TRUE)
write.csv(bg_env_table, file.path(EC.env$outputdir,
                                  paste0('background_environmental_',
                                         occur.species, '_maxent.csv')), row.names = TRUE)

## Run functions
included_layers <- enviro.data.layer[which(enviro.data.type =='continuous')]
EC_PlotContinOCC(occ_env, included_layers, outputdir=EC.env$outputdir)
EC_PlotContinBG(bg_env, included_layers, outputdir=EC.env$outputdir)

included_layers <- enviro.data.layer[which(enviro.data.type !='continuous')]
EC_PlotCategOCC(occ_env, included_layers, outputdir=EC.env$outputdir)
EC_PlotCategBG(bg_env, included_layers, outputdir=EC.env$outputdir)

#===================================================================
## MaxEnt modeling

# The Maxent function in the 'biomod2' package uses the following arguments:
# x: Predictors, this is the RasterStack with the predictor layers,
# p: Occurrence data,
# a: Background points, 
# nbg: Number of background points. Sampled randomly from the cells that are not 
#      NA in the first predictor variable. 
# This argument is ignored if background points are specified with argument a
# factors: which predictor variables, if any, should be considered as categorical?
# args: additonal arguments

# Set all Maxent arguments, based on the parameters in your main file
maxent_args <- c(
  # remove duplicate presence records. If environmental data are in grids, 
  # duplicates are records in the same grid cell or with identical coordinates.
  paste0('removeduplicates=', EC.params$removeduplicates), 
  # if the number of background points/grid cells is larger than this number, 
  # then this number of cells is chosen randomly for background points.
  paste0('maximumbackground=', EC.params$maximumbackground), 
  # add to the background any sample for which has a combination of environmental
  # values that isn't already present in the background
  paste0('addsamplestobackground=', EC.params$addsamplestobackground), 
  # add all samples to the background, even if they have combinations of 
  # environmental values that are already present in the background
  paste0('addallsamplestobackground=', EC.params$addallsamplestobackground), 
  
  # missing data
  # on model training, allow samples that have nodata values for one or more environmental variables
  paste0('allowpartialdata=', EC.params$allowpartialdata),
  
  # variable importance
  # default=FALSE; measure importance of each environmental variable by training
  # with each environmental variable first omitted, then used in isolation.
  paste0('jackknife=', EC.params$jackknife),
  
  # random seed
  # if selected, a different random seed will be used for each run, so a different
  # random test/train partition will be made and a different subset of the background
  # will be used, if applicable.
  paste0('randomseed=', EC.params$randomseed),
  
  # prevalence
  # default: probability of presence at ordinary occurrence points. 
  # See Elith et al. Diversity and Distributions, 2011 for details
  paste0('defaultprevalence=', EC.params$defaultprevalence),
  
  # train/test settings
  # percentage of presence localities to be randomly set as test points
  paste0('randomtestpoints=', EC.params$randomtestpoints), 
  # number of rep runs to do when cross-validating/bootstrapping/sampling with replacement runs.
  paste0('replicates=', EC.params$replicates), 
  # if replicates > 1, do multiple runs of this type 
  paste0('replicatetype=', EC.params$replicatetype), 
  # stop training after this many iterations of the optimization algorithm
  paste0('maximumiterations=', EC.params$maximumiterations),
  # stop training when the drop in log loss per iteration drops below this number 
  paste0('convergencethreshold=', EC.params$convergencethreshold),
  
  # feature selection
  # automatically select features classes to use, based on number of training samples
  paste0('autofeature=', EC.params$autofeature), 
  paste0('linear=', EC.params$linear),  # allow linear features to be used
  paste0('quadratic=', EC.params$quadratic),  # allow quadratic features to be used
  paste0('product=', EC.params$product),  # allow product features to be used
  paste0('threshold=', EC.params$threshold),  # allow threshold features to be used
  paste0('hinge=', EC.params$hinge),  # allow hinge features to be used
  
  # feature settings
  # number of samples at which product and threshold features start being used
  paste0('lq2lqptthreshold=', EC.params$lq2lqptthreshold),
  # number of samples at which quadratic features start being used
  paste0('l2lqthreshold=', EC.params$l2lqthreshold),
  # number of samples at which hinge features start being used
  paste0('hingethreshold=', EC.params$hingethreshold),
  
  # regularization settings
  # multiply automatic regularization parameters by this number
  # a higher number gives a more spread-out distribution.
  paste0('betamultiplier=', EC.params$betamultiplier),
  # parameter to be applied to all threshold features; negative value enables automatic setting
  paste0('beta_threshold=', EC.params$beta_threshold),
  # to be applied to all categorical features; negative value enables automatic setting
  paste0('beta_categorical=', EC.params$beta_categorical), 
  # to be applied to all linear, quadratic and product features; negative value enables automatic setting
  paste0('beta_lqp=', EC.params$beta_lqp),
  #to be applied to all hinge features; negative value enables automatic setting
  paste0('beta_hinge=', EC.params$beta_hinge),
  
  # outputs - These are not shown in the UI, so unable to be changed by user
  'responsecurves=TRUE',  # default=FALSE shows predicted relative probability of each environmental variable
  'responsecurvesexponent=FALSE',  # exponent (instead of logistic value) for the y axis in response curves
  'pictures=TRUE',  # create a .png image for each output grid
  'outputformat=raw',  # probabilities used in writing output grids, see Help for details
  'writeclampgrid=TRUE',  # grid with the absolute difference between prediction values with and without clamping
  # create a multidimensional environmental similarity surface (MESS)
  'writemess=TRUE',  # shows where novel climate conditions exist in the projection layers, 
                     # with the degree of novelness and the variable that is most out of range at each point
  'writeplotdata=FALSE',  # output files with the data used to make response curves
  'outputgrids=TRUE',  # if off when doing replicate runs, causes only the summary grids to be written
  'plots=TRUE',  # write various plots for inclusion in .html output
  'logfile=maxent.log',  # to be used for writing debugging information about a run in output directory
 #'applythresholdrule=Fixed cumulative value 1',  # threshold rule, generating a binary outputgrid in addition 
                    #'#to the regular prediction grid. Use the full name of the threshold rule in Maxent's html output as the argument.
  'logscale=TRUE',  # all pictures of models will use a logarithmic scale for color-coding
  'writebackgroundpredictions=FALSE',  # write.csv file with predictions at background points
  'fadebyclamping=FALSE',  # reduce prediction by the difference between clamped and non-clamped output at that point

  #projection settings NB. These are not shown in the UI, so unable to be changed by user
  'extrapolate=TRUE',  # predict to regions of environmental space outside the limits encountered during training
  'doclamp=TRUE'  # apply clamping when projecting
  
)

#===================================================================
## Run model for the raw output

writeLines('Running model for raw output ...')
model_raw <- dismo::maxent(current.climate.scenario, occ_coords, bg_biased, 
                    path=file.path(EC.env$outputdir, "raw"), 
    args=maxent_args)


# Get entropy value used to transform the raw output to logistc and cloglog
lambdas <- rmaxent::parse_lambdas(model_raw)

# Project the raw model to model-fitting grids
pred.raw <- rmaxent::predict(model_raw, current.climate.scenario, args=c(outputformat="raw"))

# Transform the raw output to logistic/cloglog output, used for the habitat suitability maps
pred.cloglog2 <- rmaxent::to_cloglog(pred.raw, from="raw", lambdas$entropy)
pred.logistic2 <- rmaxent::to_logistic(pred.raw, from="raw", lambdas$entropy)

#===================================================================
## Maxent model results

# Habitat suitability map for each output format
fname = paste0("proj_current_", occur.species, "_maxent_cloglog")
writeRaster(pred.cloglog2, filename=file.path(EC.env$outputdir, fname), 
            format = "GTiff", overwrite=TRUE)
png(filename=file.path(EC.env$outputdir, paste0(fname, ".png")), width=800, height=800)
title = "Habitat Suitability projection (cloglog)"
plot(pred.cloglog2, main=title, xlab='longitude', ylab='latitude')
dev.off()

fname = paste0("proj_current_", occur.species, "_maxent_logistic")
writeRaster(pred.logistic2, filename=file.path(EC.env$outputdir, fname),
            format = "GTiff", overwrite=TRUE)
png(filename=file.path(EC.env$outputdir, paste0(fname, ".png")), width=800, height=800)
title = "Habitat Suitability projection (logistic)"
plot(pred.logistic2, main=title, xlab='longitude', ylab='latitude')
dev.off()


# Response curves

# Generate response curves for cloglog and logistic output
# This script converts the output maps and figures to log and clog log (raw >log> cloglog)
# These are the maxEnt models for the logistic and cloglog outputs
writeLines('Running model for cloglog output ...')
model_rcCloglog <- maxent(current.climate.scenario, occ_coords, bg_biased, 
                          path=file.path(EC.env$outputdir, "cloglog"), 
                          args=maxent_args, outputformat="cloglog")
writeLines('Running model for logistic output ...')
model_rcLogistic <- maxent(current.climate.scenario, occ_coords, bg_biased, 
                           path=file.path(EC.env$outputdir, "logistic"), 
                          args=maxent_args, outputformat="logistic")

# Copy the ROC file and omission file from cloglog results
tfile = file.path(EC.env$outputdir, "raw", "plots", "species_roc.png")
if (file.exists(tfile)) {
  file.copy(tfile, file.path(EC.env$outputdir, paste0("roc_", occur.species, "_maxent.png")))
}

tfile = file.path(EC.env$outputdir, "raw", "plots", "species_omission.png")
if (file.exists(tfile)) {
  file.copy(tfile, file.path(EC.env$outputdir, paste0("omcom_", occur.species, "_maxent.png")))
}

# Plot response curves
EC_PlotResponseCurve(model_rcCloglog, names(current.climate.scenario), occur.species, "cloglog")
EC_PlotResponseCurve(model_rcLogistic, names(current.climate.scenario), occur.species, "logistic")

# Create a Maxent results table
maxentResults <- read.csv(file.path(EC.env$outputdir, "raw", "maxentResults.csv"))
maxentResults <- t(maxentResults)
head(maxentResults)
write.csv(maxentResults, file.path(EC.env$outputdir, paste0("maxent_results_", occur.species, ".csv")))

# Multivariate Environmental Similarity Surface (MESS)
# NB. Default is full=TRUE however this returns a raster brick with all variables
# Changed to full=FALSE to get only MESS result
mess <- suppressWarnings(mess(current.climate.scenario, occ_env, full=FALSE))
writeRaster(x=mess, file.path(EC.env$outputdir, paste0("mess_", occur.species, "_maxent.tif")),
            format="GTiff", overwrite=TRUE)
png(filename=file.path(EC.env$outputdir, paste0("mess_", occur.species, "_maxent.png")), 
    width=800, height=800)
plot(mess, col=terrain.colors(10, alpha = 1), xlab='longitude', ylab='latitude')
dev.off()

## Limiting Factors Map which is categorical data
limiting_factors <- limiting(current.climate.scenario, model_raw)
png(filename=file.path(EC.env$outputdir, paste0("limitfactors_", occur.species, "_maxent.png")),
    width=800, height=800)
lf_plot = levelplot(limiting_factors)
print(lf_plot)
dev.off()

#===================================================================
## Maxent model evaluation
## Random NULL model

# Usage:
# x = data frame with environmental predictor values for occurrences
# model = model function that creates a model of class 'DistModel', example uses bioclim
# n = sample size
# rep = number of repetitions
# pa = boolean (TRUE or FALSE), presence-only or presence-background

writeLines('Random NULL model ...')
nr <- nullRandom(as.data.frame(occ_env), maxent, n=min(25, nrow(occ_env)-1), 
                 rep=EC.params$randomnullrep, pa=TRUE)

# Determine the 95% confidence interval (CI)
CI <- nr[[round(0.95*EC.params$randomnullrep)]]
print(CI)

# save the model AUC and null model AUC
auc_df = data.frame('Training Model AUC'=maxentResults['Training.AUC',], 'NULL MODEL AUC'=CI@auc)
write.csv(auc_df, file.path(EC.env$outputdir, paste0("nullmodel_", occur.species, "_maxent.csv")))

#===================================================================
##  Evaluation statistics

equal_sens_spec_threshold <- EC_SummaryStats(occ_coords,
                                                      bg_biased, model_raw, 
                                                      current.climate.scenario, occur.species)
writeLines(paste0('equal_sens_spec_threshold = ', equal_sens_spec_threshold))

# Area/Extent of Occuppancy (AOO/EOO)
#EC_AreaExtentOccupancy(occur, occur.species)

# generate species threshold map based on cloglog projection 
proj_file <- file.path(EC.env$outputdir, paste0("proj_current_",
                                                occur.species, "_maxent_cloglog.tif"))
proj_raster <- EC.raster.load(proj_file)
EC_ThresholdMap(proj_raster, equal_sens_spec_threshold)

#===================================================================
## Projection over current climate scenario without constraint 
#  NOTE: only if all env data layers are continuous

if (enviro.data.genUnconstraintMap && 
   all(enviro.data.type == 'continuous') && 
   (!is.null(enviro.data.constraints) || enviro.data.generateCHull)) {

   writeLines("Doing projection over unconstraint climate scenario ...")
   pred.raw.unconstraint <- predict(model_raw, current.climate.scenario.orig, 
                                    args=c(outputformat="raw"))

   # raw output to logistic/cloglog output;used to produce the habitat suitability maps
   pred.cloglog2_unconstraint <- rmaxent::to_cloglog(pred.raw.unconstraint, 
                                                     from="raw", lambdas$entropy)
   pred.logistic2_unconstraint <- rmaxent::to_logistic(pred.raw.unconstraint,
                                                     from="raw", lambdas$entropy)

  # habitat suitability map for each output format
  fname = paste0("proj_current_", occur.species, "_maxent_cloglog_unconstraint")
  writeRaster(pred.cloglog2_unconstraint, filename=file.path(EC.env$outputdir, 
                                                             fname), 
              format = "GTiff", overwrite=TRUE)
  fname = paste0("proj_current_", occur.species, "_maxent_logistic_unconstraint")
  writeRaster(pred.logistic2_unconstraint, filename=file.path(EC.env$outputdir, 
                                                              fname),
              format = "GTiff", overwrite=TRUE)
}

writeLines('Done')