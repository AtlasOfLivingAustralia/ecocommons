############################################################
###       EcoCommons script to run SDM algorithms        ###
############################################################
##
## Author details: EcoCommons Platform 
## Contact details: comms@ecocommons.org.au
## Copyright statement: This script is the product of EcoCommons platform. Please
## run "Citation" to check how to cite this R script, etc, etc
## Date: September 2021
##
## Script and data info:
## In this script you will: 
## - Set up parameters that will allow you to run multiple algorithms for 
##    species distribution modelling (SDMs)
##
## - Run the SDMs with selected algorithms

# You need devtools to install the ecocommons package from GitHub
install.packages ("devtools")
library (devtools)
devtools::install_github ("AtlasOfLivingAustralia/ecocommons")
library(ecocommons)


# 1. Load and edit your dataset

# The 'param.json' file used here is generated on the EcoCommons platform
source_file <- EC_read_json (file="~/Documents/GitHub/ecocommons_ALA/inst/variables/test_constraint_map_ann.json")

print.model_parameters (source_file)

EC.params <- source_file$params  # data file with species, environment data and parameters

EC.env <- source_file$env  # set workplace environment


# Print out parameters used
parameter.print(EC.params)

# 2. Set working directory (script runner takes care of it)
setwd(EC.env$workdir)

# 3. Provide simple names for all the parameters
  # parameter.as.spring
  # parameter.print

# 4. Set data for modelling
response_info <- EC_build_response (EC.params)  # species data

predictor_info <- EC_build_predictor (EC.params)  # environmental data

constraint_info <- EC_build_constraint (EC.params)  # constraint area data

dataset_info <- EC_build_dataset (predictor_info, constraint_info,
                                  response_info)  # read data and constraint to region of interest

# Compute the model using biomod2 format; create pseudo-absence if no true absence is given
model_compute <- EC_format_biomod2 (true.absen               = dataset_info$absen,
                                    pseudo.absen.points      = dataset_info$pa_number_point,
                                    pseudo.absen.strategy    = 'random',
                                    pseudo.absen.disk.min    = EC.params$pa_disk_min,
                                    pseudo.absen.disk.max    = EC.params$pa_disk_max,
                                    pseudo.absen.sre.quant   = EC.params$pa_sre_quant,
                                    climate.data             = dataset_info$current_climate,
                                    occur                    = dataset_info$occur,
                                    species.name             = response_info$occur_species,
                                    save.pseudo.absen        = TRUE,
                                    save.env.absen           = TRUE,
                                    save.env.occur           = TRUE,
                                    generate.background.data = FALSE)

# OPTION 2: return a json here; than can implement each algorithm separately


#OPTION 1: If it is just one JSON, read below:

if ( 'ANN' | 'CTA' | 'GBM' | 'RF' | 'FDA' | 'GAM' | 'GLM' | 'MARS')
  # machine learning and statistical regression
  
  # specific parameters restricted to each algorithm
  model_options_ann <- EC_options_ann (EC.params$ann) # if all algorithms are in the same list

  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (ANN = model_options$ann,
                                                    GLM = model_options$glm)  # if a list is created


  # Use define model options (created on model_compute) and compute the model
  # uses biomod2 model options
  model_sdm <-
    biomod2::BIOMOD_Modeling (data               = model_compute,
                              models             = model_algorithm,  # can be a list of multiple algorithms
                              models.options     = model_options,
                              NbRunEval          = model_options_algorithm$NbRunEval,
                              DataSplit          = model_options_algorithm$DataSplit,
                              Yweights           = model_options_algorithm$Yweights,
                              Prevalence         = model_options_algorithm$Prevalence,
                              VarImport          = model_options_algorithm$VarImport,
                              models.eval.meth   = model_options_algorithm$biomod_eval_meth,
                              SaveObj            = TRUE,
                              rescal.all.models  = model_options_algorithm$rescal_all_models,
                              do.full.models     = model_options_algorithm$do_full_models,
                              modeling.id        = model_options_algorithm$model_id)
  
  # save out the model object
  EC_save (model_sdm, name = "model.object.RData")
  




#########____________EXAMPLES______________#########
  
# GEOGRAPHICAL TYPE ALGORITHMS
koala_circles <- EC_modelling_circles (EC.params,response_info,
                                       predictor_info, dataset_info)

koala_convex_hull <- EC_modelling_convhull (EC.params,response_info,
                                            predictor_info, dataset_info)

koala_geo_distance <- EC_modelling_geodist (EC.params,response_info,
                                            predictor_info, dataset_info)

koala_inversed_distance <- EC_modelling_geoidw (EC.params,response_info,
                                                predictor_info, dataset_info)

koala_voronoi_hull <- EC_modelling_voronoihull (EC.params,response_info,
                                                predictor_info, dataset_info)


# STATISTICAL REGRESSION ALGORITHMS
# uses EC_utilities_algorithms + individual functions
koala_FDA <- EC_modelling_fda (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GAM <- EC_modelling_gam (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GLM <- EC_modelling_glm (EC.params,response_info,
                               predictor_info, dataset_info)

koala_mars <- EC_modelling_mars (EC.params,response_info,
                                 predictor_info, dataset_info)



# MACHINE LEARNING TYPE ALGORITHMS
# uses EC_utilities_algorithms + individual functions
koala_ANN <- EC_modelling_ann (EC.params,       # EC.params
                               response_info,   # from EC_build_response
                               model_compute)   # from EC_build_dataset 

koala_BRT <- EC_modelling_brt (EC.params,response_info,
                               predictor_info, dataset_info)

koala_CTA <- EC_modelling_cta (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GBM <- EC_modelling_gbm (EC.params,response_info,
                               predictor_info, dataset_info)

koala_random_forest <- EC_modelling_rf (EC.params,response_info,
                                        predictor_info, dataset_info)


# MAXENT
# uses EC_utilities_maxent + individual functions
koala_maxent <- EC_modelling_maxent (EC.params,response_info,
                                     predictor_info, dataset_info)

# PROFILE TYPE ALGORITHMS
# uses EC_utilities_algorithms + individual functions
koala_bioclim <- EC_modelling_bioclim (EC.params,response_info,
                                       predictor_info, dataset_info)

koala_SRE <- EC_modelling_sre (EC.params,response_info,
                               predictor_info, dataset_info)



## DONE!