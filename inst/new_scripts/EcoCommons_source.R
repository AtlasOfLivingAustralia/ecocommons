############################################################
###       EcoCommons script to run SDM algorithms        ###
############################################################
##
## Author details: EcoCommons Platform 
## Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons platform. Please
## run "Citation" to check how to cite this R script, etc, etc
## Date : ??
##
## Script and data info:
## In this script you will 
## 1) Set up parameters that will allow you to run multiple algorithms for 
##    species distribution modelling (SDMs);
## 2) Run the SDMs with selected algorithms
##

# Load and edit your dataset

# The 'param.json' file used here is generated on EcoCommons python script
source_file <- EC_read_json(file="~/Documents/GitHub/ecocommons_ALA/inst/variables/test_constraint_map.json")

print.model_parameters(source_file)

EC.params <- source_file$params  # data file with species, environment data and parameters

EC.env <- source_file$env  # set worksplace environment


# Set working directory (script runner takes care of it)
setwd(EC.env$workdir)

# Set a temporary directory to raster - we do this because raster sometimes makes
# temp files (e.g. when cropping).Might want to make this configurable - e.g. we
# might want to control maxmemory, and/or other raster options
# rasterOptions(tmpdir=paste(EC.env$workdir,"raster_tmp",sep="/"))


# Set data for modelling
response_info <- EC_build_response(EC.params)  # species data

predictor_info <- EC_build_predictor(EC.params)  # environmental data

constraint_info <- EC_build_constraint(EC.params)  # constraint area data

dataset_info <- EC_build_dataset(response_info, predictor_info,
                              constraint_info)

# Choose algorithm to run
# Include a list of available algorithms to run

#===================================================================
# EXAMPLES - not for the user, just for testing now

# GEOGRAPHICAL TYPE ALGORITHMS
# uses EC_utilities_geographical + individual functions
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

############
## STOP HERE! 
############

# STATISTICAL REGRESSION ALGORITHMS
# uses EC_utilities_prof+stat_regression + individual functions
koala_FDA <- EC_modelling_fda (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GAM <- EC_modelling_gam (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GLM <- EC_modelling_glm (EC.params,response_info,
                               predictor_info, dataset_info)

koala_mars <- EC_modelling_mars (EC.params,response_info,
                                 predictor_info, dataset_info)



# PROFILE TYPE ALGORITHMS
# uses EC_utilities_prof+stat_regression + individual functions
koala_bioclim <- EC_modelling_bioclim (EC.params,response_info,
                                       predictor_info, dataset_info)

koala_SRE <- EC_modelling_sre (EC.params,response_info,
                               predictor_info, dataset_info)



# MACHINE LEARNING TYPE ALGORITHMS
# uses EC_utilities_machine_learning + individual functions
koala_ANN <- EC_modelling_ann (EC.params,response_info,
                               predictor_info, dataset_info)

koala_BRT <- EC_modelling_brt (EC.params,response_info,
                               predictor_info, dataset_info)

koala_CTA <- EC_modelling_cta (EC.params,response_info,
                               predictor_info, dataset_info)

koala_GBM <- EC_modelling_gbm (EC.params,response_info,
                               predictor_info, dataset_info)

koala_maxent <- EC_modelling_maxent (EC.params,response_info,
                                     predictor_info, dataset_info)

koala_random_forest <- EC_modelling_rf (EC.params,response_info,
                                        predictor_info, dataset_info)



#===================================================================
## Species trait data
## still need to work on it

# Link to input dataset csv file
trait.data.filename <- EC.params$traits_dataset$filename
# Link to variable names of input dataset
trait.data.params <- EC.params$traits_dataset_params
# Read in the trait data
trait.data <- read.csv(trait.data.filename)
# Get the species
trait.species <-EC.params$species


writeLines("Done!")
