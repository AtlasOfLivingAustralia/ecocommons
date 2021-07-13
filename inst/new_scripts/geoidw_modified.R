################################################################################
###   EcoCommons script to Inverse Distance Weighted Model - using biomod2   ###
################################################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs Inverse Distance Weighted Model to predicts a species 
##  probabilities in unknown locations by averaging values of nearby occurrences.
##  It is a geographical method for Species Distribution Modelling. It uses 
##  the biomod2 package on R to set up parameters, dismo functions, and the 
##  input data generated on the EcoCommons platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## Load and edit your dataset with parameters specific to Inverse Distance Weighted Model



#additional parameters for projecting circles
opt.tails <- EC.params$tails # default "both"; use to ignore the left or right tail of the percentile distribution ("both", "low", "high"
opt.ext <- NULL #an optional extent object to limit the prediction to a sub-region of 'x'
projection.name <- "current"
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_geoIDW", occur.species), 
                          sprintf("%s_geoIDW_%s", occur.species, EC.params$subset))


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth <- c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY",
                             "BIAS", "POD", "CSI", "ETS")
# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy <- c(dismo.eval.method, biomod.models.eval.meth)
# TODO: these functions are used to evaluate the model ... configurable?


# Determine the number of pseudo absence points from pa_ratio
pa_ratio <- EC.params$pa_ratio
pa_number_point <- 0
if (pa_ratio > 0) {
  pa_number_point <- floor(pa_ratio * nrow(occur))
}

# Format the data as in biomod2. This will also generate the psedo absence points.
biomod2.data <- EC_FormatDataBIOMOD2(true.absen       = absen,
                                  pseudo.absen.points    = pa_number_point,
                                  pseudo.absen.strategy  = EC.params$pa_strategy,
                                  pseudo.absen.disk.min  = EC.params$pa_disk_min,
                                  pseudo.absen.disk.max  = EC.params$pa_disk_max,
                                  pseudo.absen.sre.quant = EC.params$pa_sre_quant,
                                  climate.data           = current.climate.scenario,
                                  occur                  = occur,
                                  species.name           = occur.species,
                                  save.env.absen         = FALSE,
                                  save.env.occur         = FALSE,
                                  species_algo_str       = species_algo_str)

#===================================================================
## GEOIDW - Inverse Distance Weighted interpolation

# geoIDW(p, a, ...)
# p presence points; two column matrix, data.frame or SpatialPoints* object
# a absence points; must be of the same class as 'p'
# ... none implemented

model.sdm <- geoIDW(p=occur, a=absen)
# save out the model object
EC_Save(model.sdm, EC_OutfileName(filename="model.object",
                                  id_str=species_algo_str, ext="RData"))

# predict for given climate scenario
model.proj <- predict(model.sdm, current.climate.scenario@layers[[1]], mask=TRUE)


# remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# save output
EC_SaveModelProjection(model.proj, projection.name, occur.species, species_algo_str)
