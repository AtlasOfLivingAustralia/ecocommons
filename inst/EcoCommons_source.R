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
devtools::install_github ("AtlasOfLivingAustralia/ecocommons") # CHECK NOW PSEUDO-ABSENCE BRANCH
library(ecocommons)


# 1. Load and edit your dataset

# The 'param.json' file used here is generated on the EcoCommons platform
source_file <- EC_read_json (file="~/Documents/GitHub/ecocommons/inst/variables/test_constraint_map_ann.json")

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

## check access to build_functions 
response_info <- EC_build_response (EC.params)  # species data

predictor_info <- EC_build_predictor (EC.params)  # environmental data

constraint_info <- EC_build_constraint (EC.params)  # constraint area data

dataset_info <- EC_build_dataset (predictor_info, constraint_info,
                                  response_info)  # read data and constraint to region of interest

# Compute the model using biomod2 format; create pseudo-absence if no true absence is given

#think about a better name
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



#########____________EXAMPLES______________#########
# For testing, not for users
  
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