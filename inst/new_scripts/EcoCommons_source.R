############################################################
###        EcoCommons script to define parameters        ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
## In this script you will set up parameters that will allow you to run
## multiple algorithms for species distribution modelling (SDMs). To start, you
## will need 1) install and load all necessary packages, 2) load and edit
## your dataset


## Load and edit your dataset

# The 'param.json' file used here is generated on EcoCommons python script
source_file <- EC_ReadJson(file="~/Documents/GitHub/ecocommons_ALA/inst/variables/test_script_circles2.json")
print.model_parameters(source_file)
a<-source_file$params  # data file with species, enviroment data and parameters
EC.env<-source_file$env  # set workspace environment


# Set working directory (script runner takes care of it)
setwd(EC.env$outputdir)

# Set a temporary directory to raster - we do this because raster sometimes makes
# temp files (e.g. when cropping).Might want to make this configurable - e.g. we
# might want to control maxmemory, and/or other raster options
# rasterOptions(tmpdir=paste(EC.env$workdir,"raster_tmp",sep="/"))

# Use seed supplied, if any. Otherwise generate a random seed.
seed <- a$random_seed
if (is.null(seed)) {
  seed = runif(1, -2^31, 2^31-1)
}
seed <- as.integer(seed)
set.seed(seed)
a["random_seed"] <- seed

# Set data for modelling
response_info <- build_response(a)  # species data
predictor_info <- build_predictor(a,response_info)  # environmental data



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
