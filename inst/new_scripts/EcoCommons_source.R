############################################################
###        EcoCommons script to define parameters        ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of BCCVL etc.
## Date : ??
## Script and data info:
## In this script you will set up parameters that will allow you to run
## multiple algorithms for species distribution modelling (SDMs). To start, you
## will need 1) install and load all necessary packages, 2) load and edit
## your dataset


# First, setup a R environment on host computer, defining a CRAN repository
# You can jump this step if running remotely at home
r <- getOption ("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos = r)
# print warnings immediately
options(warn = 1)


# Check the required packages for distribution models
# 1. print out list of installed packages
write.table(installed.packages() [,c("Package", "Version", "Priority")],
            row.names=FALSE)

# 2. check if libraries are installed, install if necessary and then load them
necessary <- c("biomod2","ConR","dismo", "dplyr", "gdalUtils", "gridExtra",
               "raster","rasterVis","rgeos","rgdal", "rjson","sp")

### packages ‘ConR’, ‘dev2’ are not available for this version of R

installed <- necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >= 1) {
  install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
  library(lib,character.only = T) #load the libraries
}

#===================================================================
## Load and edit your dataset

# The 'param.json' file used here is generated on EcoCommons python script
params <- rjson::fromJSON(file="~/Documents/BCCVL_scripts/test_script.json")
EC.params <- params$params
EC.env <- params$env
rm(params)

# Set working directory (script runner takes care of it)
setwd(EC.env$outputdir)
# Set a temporary directory to raster - we do this because raster sometimes makes
# temp files (e.g. when cropping).Might want to make this configurable - e.g. we
# might want to control maxmemory, and/or other raster options
rasterOptions(tmpdir=paste(EC.env$workdir,"raster_tmp",sep="/"))

# Use seed supplied, if any. Otherwise generate a random seed.
seed <- EC.params$random_seed
if (is.null(seed)) {
  seed = runif(1, -2^31, 2^31-1)
}
seed <- as.integer(seed)
set.seed(seed)
EC.params["random_seed"] <- seed

# The lon/lat of the observation records -- minimum of 2 column matrix
occur.data <- EC.params$species_occurrence_dataset$filename
occur.species <- EC.params$species_occurrence_dataset$species
month.filter <- EC.params$species_filter #if user has any filters set for the analyses

# The lon/lat of the background / absence points -- 2 column matrix
absen.data <- EC.params$species_absence_dataset$filename

# OPTIONAL: bias file from user
bias.data <- EC.params$bias_dataset$filename
if (is.null(bias.data) || bias.data == '') {
  bias.data = NULL
} else {
 bias.data = readLines(EC.params$bias_dataset$filename)
}


# The current environment data -- list of layers' file names
enviro.data.current <- lapply(EC.params$environmental_datasets,
                              function(x) {
                                fname = x$filename
                                return(fname)
                              }
)

# Type name in terms of continuous or categorical
enviro.data.type <- lapply(EC.params$environmental_datasets, function(x) x$type)
# Layer names for the current environmental layers used
enviro.data.layer <- lapply(EC.params$environmental_datasets, function(x) x$layer)

# OPTIONAL- you can add a geographic constraints
enviro.data.constraints <- EC.params$modelling_region
if (is.null(enviro.data.constraints) || enviro.data.constraints == '') {
  enviro.data.constraints = NULL
} else {
  enviro.data.constraints = readLines(EC.params$modelling_region$filename)
}

# Indicate if you wish to generate and apply convex-hull polygon of occurrence
# dataset to constraint
enviro.data.generateCHull <- ifelse(is.null(EC.params$generate_convexhull),
                                    FALSE, as.logical(EC.params$generate_convexhull))

# Indicate whether to generate unconstrained map or not. TRUE by default
enviro.data.genUnconstraintMap <- ifelse(is.null(EC.params$unconstraint_map),
                                         TRUE, as.logical(EC.params$unconstraint_map))
# Resampling (up / down scaling) if scale_down is TRUE, return 'lowest'
enviro.data.resampling <- ifelse(is.null(EC.params$scale_down) ||
                                   as.logical(EC.params$scale_down),
                                 'highest', 'lowest')

# Read current climate data
current.climate.scenario <- EC_EnviroStack(enviro.data.current,
                                           enviro.data.type, enviro.data.layer,
                                           resamplingflag=enviro.data.resampling)

# Read in the necessary observation, background and environmental data
occur <- EC_ReadSp(occur.data, month.filter) #read in the observation lon/lat
absen <- EC_ReadSp(absen.data, month.filter) #read in the absence lon/lat

# Geographically constrained modelling
if (!is.null(enviro.data.constraints) || enviro.data.generateCHull) {
  constrainedResults <- EC_SDMGeoConstrained(current.climate.scenario, occur,
                                             absen, enviro.data.constraints,
                                             enviro.data.generateCHull);

  # Save a copy of the climate dataset
  current.climate.scenario.orig <- current.climate.scenario
  current.climate.scenario <- constrainedResults$raster
  occur <- constrainedResults$occur
  absen <- constrainedResults$absen
}


#===================================================================
## Species trait data


# Link to input dataset csv file
trait.data.filename = EC.params$traits_dataset$filename
# Link to variable names of input dataset
trait.data.params = EC.params$traits_dataset_params
# Read in the trait data
trait.data = read.csv(trait.data.filename)
# Get the species
trait.species =EC.params$species


writeLines("Done!")
