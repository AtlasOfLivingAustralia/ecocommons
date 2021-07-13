############################################################
###      EcoCommons script to circles - using Dismo      ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs Circles to predict a species within a set radius 
##  around occurrence records. It is a geographical method for Species 
##  Distribution Modelling. It uses the biomod2 package on R to set up 
##  parameters, dismo functions, and the input data generated on the  
##  EcoCommons platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## Load and edit your dataset with parameters specific to Circles

opt.d <- EC.params$d # radius around circles; if not specified, computed from mean inter-point distance 

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
## Circles Modelling

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
EC_Save(model.sdm, EC_OutfileName(filename="model.object",
                                            id_str=species_algo_str, ext="RData"))

# Predict for given climate scenario
model.proj <- predict(model.sdm, current.climate.scenario@layers[[1]], mask=TRUE)

# Remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# Save output
EC_SaveModelPrediction(model.proj, projection.name, occur.species, species_algo_str)
