############################################################
###      EcoCommons script to circles - using Dismo      ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs the circles algorithm (ref)
##  gegraphical; for Species Distribution Modelling, using the 
##  biomod2 package on R and the input data generated on the EcoCommons 
##  platform.

#===================================================================
## Load and edit your dataset

# The lon/lat of the observation records -- minimum of 2 column matrix
occur.data <- EC.params$species_occurrence_dataset$filename
occur.species <- EC.params$species_occurrence_dataset$species
month.filter <- EC.params$species_filter #if user has any filters set for the analyses

# The lon/lat of the background / absence points -- 2 column matrix
absen.data <- EC.params$species_absence_dataset$filename

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


# Additional parameters for 'Circles' algorithm
opt.d <- EC.params$d # radius around circles, if not specified,computed from mean inter-point distance 

# Use to ignore the left or right tail of the percentile distribution ("both", "low", "high")
opt.tails <- EC.params$tails # default "both"
opt.ext <- NULL # an optional extent object to limit the prediction to a sub-region of 'x'
projection.name <- "current"
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_circle", occur.species), 
                          sprintf("%s_circle_%s", occur.species, EC.params$subset))

# Model accuracy statistics
# From dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth <- c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY", 
                             "BIAS", "POD", "CSI", "ETS")
# Model accuracy statistics; combine stats from dismo & biomod2 for consistent output
model.accuracy <- c(dismo.eval.method, biomod.models.eval.meth)
# TODO: these functions are used to evaluate the model ... configurable?


# Read current climate data
current.climate.scenario <- EC_EnviroStack(enviro.data.current, 
                                           enviro.data.type, enviro.data.layer, 
                                           resamplingflag=enviro.data.resampling)

# Read in the necessary observation, background and environmental data
occur <- EC_SpRead(occur.data, month.filter) #read in the observation lon/lat
absen <- EC_SpRead(absen.data, month.filter) #read in the absence lon/lat

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


# Circles algorithm does not need pseudo-absence points
pa_number_point <- 0

# format the data as in biomod2. This will also generate the psedo absence points.
biomod2.data <- EC_FormatDataBIOMOD2(true.absen       = absen,
                                  pseudo.absen.points    = pa_number_point,
                                  pseudo.absen.strategy  = EC.params$pa_strategy,
                                  pseudo.absen.disk.min  = EC.params$pa_disk_min,
                                  pseudo.absen.disk.max  = EC.params$pa_disk_max,
                                  pseudo.absen.sre.quant = EC.params$pa_sre_quant,
                                  climate.data           = current.climate.scenario,
                                  occur                  = occur,
                                  species.name           = occur.species,
                                  save.pseudo.absen      = FALSE,
                                  save.env.occur         = FALSE,
                                  species_algo_str       = species_algo_str)


#===================================================================
# Circles Modelling

# circles(p, d, ...)
# p point locations (presence), two column matrix, data.frame or 
# SpatialPoints* object d the radius of each circle in meters; a single number
# or a vector with elements corresponding to rows in 'p'; if missing the diameter
# is computed from the inter-point distance n how many vertices in the circle?
# default is 360 lonlat logical, are these longitude/latitude data? 
# r radius of the earth; relevant for longitude/latitude data; default is 6378137 m

# Run circles with occurrence data
if (is.null(opt.d)) {
  model.sdm <- circles(p=occur, lonlat=TRUE)
} else {
  model.sdm <- circles(p=occur, d=opt.d, lonlat=TRUE)
}

# Save out the model object
EC_Save(model.sdm, bccvl.format.outfilename(filename="model.object",
                                            id_str=species_algo_str, ext="RData"))

# Predict for given climate scenario
model.proj <- predict(model.sdm, current.climate.scenario@layers[[1]], mask=TRUE)

# Remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# Save output
EC_SaveModelPrediction(model.proj, projection.name, occur.species, species_algo_str)
