#' Function to run an Classification Tree Algorithm (maxent) within EcoCommons. It
#' is a machine learning model that repeatedly splits the datasets into groups
#' based on a threshold value of one of the env variables
#'
#' @param a a list created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom dismo maxent
#' @importFrom dismo nullRandom
#' @importFrom raster extract
#' @importFrom raster writeRaster
#' @importFrom rmaxent limiting
#' @importFrom rmaxent parse_lambdas
#' @importFrom rmaxent to_cloglog
#' @importFrom rmaxent to_logistic
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#'
#' @export EC_modelling_maxent


EC_modelling_maxent <- function (a,# EC.params
                                 response_info,  # from EC_build_response
                                 predictor_info,  # from EC_build_predictor
                                 dataset_info) {  # from EC_build_dataset

  # Set parameters to perform modelling
  model_algorithm <- 'maxent'

  # Determine the number of background points
  pa_number_point <- a$maximumbackground

  # Generate background data if not supplied
  if (is.null(absen)) {
    absen <- EC_create_background (predictor_info$layer[[1]], pa_number_point)
  }

  # Generate biased background data if a bias file is supplied
  bg_biased = absen
  if (!is.null(response_info$bias_data)) {
    bg_biased = EC_create_biased_background (absen,
                                             response_info$bias_data,
                                            pa_number_point)
  }

  ## Create data exploration outputs
  ## Extract predictor variable values for occurrence and background data
  occ_coords <- subset(occur[,1:2]) # need to have a dataframe with only lat and lon to extract predictor variable values
  occ_env <- raster::extract (dataset_info$current_climate, occ_coords)
  bg_env <- raster::extract (dataset_info$current_climate, bg_biased)

  # Save tables
  occ_env_table <- data.frame(cbind(occur, occ_env))
  bg_env_table <- data.frame(cbind(bg_biased, bg_env))
  write.csv(occ_env_table, file.path(EC.env$outputdir,
                                     paste0('occurrence_environmental_',
                                            response_info$occur_species,
                                            '_maxent.csv')),
            row.names = TRUE)

  write.csv(bg_env_table, file.path(EC.env$outputdir,
                                    paste0('background_environmental_',
                                           response_info$occur_species,
                                           '_maxent.csv')),
            row.names = TRUE)


  ## Run functions for continuous data
  included_layers = enviro.data.layer [which(predictor_info$type =='continuous')]
  EC_plot_contin_OCC (occ_env, included_layers, outputdir=EC.env$outputdir)
  EC_plot_contin_BG (bg_env, included_layers, outputdir=EC.env$outputdir)

  ## Run functions for categorical data
  included_layers = enviro.data.layer[which(predictor_info$type =='categorical')]
  EC_plot_categ_OCC (occ_env, included_layers, outputdir=EC.env$outputdir)
  EC_plot_categ_BG (bg_env, included_layers, outputdir=EC.env$outputdir)

  # specific parameters to run MaxEnt algorithm
  model_options_maxent <- EC_options_maxent (a)


  #--------------------------------------------------------------------------------------
  # Run model for the raw output
  writeLines('Running model for raw output ...')
  model_raw <- dismo::maxent (dataset_info$current_climate, occ_coords, bg_biased,
                              path = file.path(EC.env$outputdir, "raw"),
                              args = maxent_args)


  # Get entropy value used to transform the raw output to logistc and cloglog
  lam <- rmaxent::parse_lambdas(model_raw)

  # Project the raw model to model-fitting grids
  pred.raw <- predict(model_raw, dataset_info$current_climate, args = c(outputformat="raw"))

  # Transform the raw output to logistic/cloglog output -> this is used to produce the habitat suitability maps
  pred.cloglog2 <- rmaxent::to_cloglog (pred.raw, from="raw", lam$entropy)
  pred.logistic2 <- rmaxent::to_logistic (pred.raw, from="raw", lam$entropy)

  #------------------------------------------------------------------------------------------------------------
  ## Maxent model results
  #------------------------------------------------------------------------------------------------------------
  ## Habitat suitability map for each output format
  fname = paste0("proj_current_", response_info$occur_species, "_maxent_cloglog")
  raster::writeRaster(pred.cloglog2, filename=file.path(EC.env$outputdir, fname),
                      format = "GTiff", overwrite=TRUE)

  png(filename=file.path(EC.env$outputdir, paste0(fname, ".png")),
      width=800, height=800)
  title = "Habitat Suitability projection (cloglog)"
  plot(pred.cloglog2, main=title, xlab='longitude', ylab='latitude')
  dev.off()

  fname = paste0("proj_current_", response_info$occur_species, "_maxent_logistic")
  raster::writeRaster(pred.logistic2, filename=file.path(EC.env$outputdir, fname),
                      format = "GTiff", overwrite=TRUE)
  png(filename=file.path(EC.env$outputdir, paste0(fname, ".png")),
      width=800, height=800)
  title = "Habitat Suitability projection (logistic)"
  plot(pred.logistic2, main=title, xlab='longitude', ylab='latitude')
  dev.off()


  ## Response curves
  # Generate response curves for cloglog and logistic output
  # This script converts the output maps and figures to log and clog log (raw >log> cloglog)
  # These are the maxent models for the logistic and cloglog outputs
  writeLines('Running model for cloglog output ...')
  model_rcCloglog <- dismo::maxent (dataset_info$current_climate,
                                    occ_coords,
                                    bg_biased, path=file.path(EC.env$outputdir,
                                                              "cloglog"),
                                    args=maxent_args, outputformat="cloglog")
  writeLines('Running model for logistic output ...')

  model_rcLogistic <- maxent (dataset_info$current_climate, occ_coords,
                              bg_biased, path=file.path(EC.env$outputdir,
                                                        "logistic"),
                              args=maxent_args, outputformat="logistic")

  # Copy the roc file and omission file from cloglog results
  tfile = file.path(EC.env$outputdir, "raw", "plots", "species_roc.png")
  if (file.exists(tfile)) {
    file.copy(tfile, file.path(EC.env$outputdir, paste0("roc_",
                                                        response_info$occur_species,
                                                        "_maxent.png")))
  }

  tfile = file.path(EC.env$outputdir, "raw", "plots", "species_omission.png")
  if (file.exists(tfile)) {
    file.copy(tfile, file.path(EC.env$outputdir, paste0("omcom_",
                                                        response_info$occur_species,
                                                        "_maxent.png")))
  }

  # Plot response curves
  EC_plot_response_curve(model_rcCloglog,
                         names(dataset_info$current_climate),
                               response_info$occur_species, "cloglog")
  EC_plot_response_curve(model_rcLogistic,
                         names(dataset_info$current_climate),
                               response_info$occur_species, "logistic")

  ## Maxent results table
  maxentResults <- read.csv(file.path(EC.env$outputdir, "raw", "maxentResults.csv"))
  maxentResults <- t(maxentResults)
  head(maxentResults)
  write.csv(maxentResults, file.path(EC.env$outputdir,
                                     paste0("maxent_results_",
                                            response_info$occur_species, ".csv")))

  ## Multivariate Environmental Similarity Surface (MESS)
  # NB. Default is full=true however this returns a raster brick with all variables.
  #Changed to full=FALSE to get only MESS result.
  mess <- suppressWarnings(mess(dataset_info$current_climate, occ_env, full=FALSE))
  raster::writeRaster(x=mess, file.path(EC.env$outputdir,
                                        paste0("mess_",
                                               response_info$occur_species,
                                               "_maxent.tif")),
                      format="GTiff", overwrite=TRUE)

  png(filename = file.path(EC.env$outputdir,paste0("mess_",
                                                    response_info$occur_species,
                                                    "_maxent.png")),
      width=800, height=800)

  plot(mess, col=terrain.colors(10, alpha = 1), xlab='longitude', ylab='latitude')
  dev.off()

  ## Limiting Factors Map which is categorical data
  limiting_factors <- rmaxent::limiting (dataset_info$current_climate, model_raw)
  png(filename=file.path(EC.env$outputdir, paste0("limitfactors_",
                                                  response_info$occur_species,
                                                  "_maxent.png")),
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
  nr <- dismo::nullRandom (as.data.frame(occ_env), maxent,
                           n=min(25, nrow(occ_env)-1),
                           rep=a$randomnullrep, pa=TRUE)

  # Determine the 95% confidence interval (CI)
  CI <- nr[[round(0.95*a$randomnullrep)]]
  print(CI)

  # save the model AUC and null model AUC
  auc_df = data.frame('Training Model AUC'=maxentResults['Training.AUC',],
                      'NULL MODEL AUC'=CI@auc)
  write.csv(auc_df, file.path(EC.env$outputdir,
                              paste0("nullmodel_",
                                     response_info$occur_species,
                                     "_maxent.csv")))

  #-------------------------------------------------------------------------------
  ##  Evaluation statistics
  equal_sens_spec_threshold = EC_save_eval (occ_coords, bg_biased,
                                            model_raw,
                                            dataset_info$current_climate,
                                             response_info$occur_species)
  writeLines(paste0('equal_sens_spec_threshold = ', equal_sens_spec_threshold))

  ## Area/Extent of Occuppancy (AOO/EOO)
  #EC_area_occupancy (occur, response_info$occur_species)

  # generate species threshold map based on cloglog projection
  proj_file = file.path(EC.env$outputdir, paste0("proj_current_",
                                                 response_info$occur_species,
                                                 "_maxent_cloglog.tif"))
  proj_raster = EC_read_raster (proj_file)
  EC_threshold_map (proj_raster, equal_sens_spec_threshold)

  #--------------------------------------------------------------------------------
  # Do projection over current climate scenario without constraint only if all env data layers are continuous.
  if (enviro.data.genUnconstraintMap &&
      all(enviro.data.type == 'continuous') &&
      (!is.null(enviro.data.constraints) || enviro.data.generateCHall)) {

    writeLines("Doing projection over unconstraint climate scenario ...")
    pred.raw.unconstraint <- predict(model_raw, dataset_info$current_climate.orig, args=c(outputformat="raw"))

    # Transform the raw output to logistic/cloglog output -> this is used to produce the habitat suitability maps
    pred.cloglog2_unconstraint <- rmaxent::to_cloglog(pred.raw.unconstraint, from="raw", lam$entropy)
    pred.logistic2_unconstraint <- rmaxent::to_logistic(pred.raw.unconstraint, from="raw", lam$entropy)

    ## Habitat suitability map for each output format
    fname = paste0("proj_current_", response_info$occur_species,
                   "_maxent_cloglog_unconstraint")
    raster::writeRaster(pred.cloglog2_unconstraint,
                        filename=file.path(EC.env$outputdir, fname),
                        format = "GTiff", overwrite=TRUE)

    fname = paste0("proj_current_", response_info$occur_species,
                  "_maxent_logistic_unconstraint")
    raster::writeRaster(pred.logistic2_unconstraint,
                        filename=file.path(EC.env$outputdir, fname),
                        format = "GTiff", overwrite=TRUE)
  }

  # RETURN?? not set yet

} # end EC_modelling_maxent()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_maxent()

EC_options_maxent <- function(a){
  # Set specific parameters to run maxent algorithm
  list( # Set all Maxent arguments
    #duplicate records
    paste0('removeduplicates=', a$removeduplicates), #remove duplicate presence records. If environmental data are in grids, duplicates are records in the same grid cell, otherwise, duplicates are records with identical coordinates.

    #background records
    paste0('maximumbackground=', a$maximumbackground), #if the number of background points/grid cells is larger than this number, then this number of cells is chosen randomly for background points.
    paste0('addsamplestobackground=', a$addsamplestobackground), #add to the background any sample for which has a combination of environmental values that isn't already present in the background
    paste0('addallsamplestobackground=', a$addallsamplestobackground), #add all samples to the background, even if they have combinations of environmental values that are already present in the background

    #missing data
    paste0('allowpartialdata=', a$allowpartialdata), #during model training, allow use of samples that have nodata values for one or more environmental variables

    #variable importance
    paste0('jackknife=', a$jackknife), #default=FALSE; measure importance of each environmental variable by training with each environmental variable first omitted, then used in isolation.

    #random seed
    paste0('randomseed=', a$randomseed), #if selected, a different random seed will be used for each run, so a different random test/train partition will be made and a different random subset of the background will be used, if applicable.

    #prevalence
    paste0('defaultprevalence=', a$defaultprevalence), #default prevalence of the species: probability of presence at ordinary occurrence points. See Elith et al. Diversity and Distributions, 2011 for details

    #train/test settings
    paste0('randomtestpoints=', a$randomtestpoints), #percentage of presence localities to be randomly set aside as test points, used to compute AUC, omission, etc.
    paste0('replicates=', a$replicates), #number of replicate runs to do when cross-validating, bootstrapping or doing sampling with replacement runs.
    paste0('replicatetype=', a$replicatetype), #if replicates > 1, do multiple runs of this type:
    paste0('maximumiterations=', a$maximumiterations), #stop training after this many iterations of the optimization algorithm
    paste0('convergencethreshold=', a$convergencethreshold), #stop training when the drop in log loss per iteration drops below this number

    #feature selection
    paste0('autofeature=', a$autofeature), #automatically select which feature classes to use, based on number of training samples
    paste0('linear=', a$linear), #allow linear features to be used
    paste0('quadratic=', a$quadratic), #allow quadratic features to be used
    paste0('product=', a$product), #allow product features to be used
    paste0('threshold=', a$threshold), #allow threshold features to be used
    paste0('hinge=', a$hinge), #allow hinge features to be used

    #feature settings
    paste0('lq2lqptthreshold=', a$lq2lqptthreshold), #number of samples at which product and threshold features start being used
    paste0('l2lqthreshold=', a$l2lqthreshold), #number of samples at which quadratic features start being used
    paste0('hingethreshold=', a$hingethreshold), #number of samples at which hinge features start being used

    #regularization settings
    paste0('betamultiplier=', a$betamultiplier), #multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.
    paste0('beta_threshold=', a$beta_threshold), #regularization parameter to be applied to all threshold features; negative value enables automatic setting
    paste0('beta_categorical=', a$beta_categorical), #regularization parameter to be applied to all categorical features; negative value enables automatic setting
    paste0('beta_lqp=', a$beta_lqp), #regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting
    paste0('beta_hinge=', a$beta_hinge), #regularization parameter to be applied to all hinge features; negative value enables automatic setting

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
}

