####################################################################
###   EcoCommons script to Geographic Distance - using biomod2   ###
####################################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs Geographic Distance to predict a species by assuming 
##  that a species is more likely to be found closer to an occurrence record. 
##  It is a geographical method for Species Distribution Modelling. It uses 
## the biomod2 package on R to set up parameters, dismo functions, and the 
## input data generated on the EcoCommons platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## Load and edit your dataset with parameters specific to Geographic Distance

opt.scale <- ifelse(is.null(EC.params$scale), 1000, EC.params$scale)
opt.tails <- EC.params$tails # default "both"; use to ignore the left or right
                             # tail of the percentile distribution
opt.ext <- NULL #an optional extent object to limit the prediction to a sub-region of 'x'
projection.name <- "current"
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_geoDist", occur.species), 
                          sprintf("%s_geoDist_%s", occur.species,
                                  EC.params$subset))


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth <- c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY",
                            "BIAS", "POD", "CSI", "ETS")

# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy <- c(dismo.eval.method, biomod.models.eval.meth)

# Determine the number of pseudo absence points from pa_ratio
pa_ratio <- EC.params$pa_ratio
pa_number_point <- 0
if (pa_ratio > 0) {
  pa_number_point = floor(pa_ratio * nrow(occur))
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
                                  save.pseudo.absen      = FALSE,
                                  save.env.occur         = FALSE,
                                  species_algo_str       = species_algo_str)

# Extract occurrence and absence data
coord <- biomod2.data@coord
occur <- coord[c(which(biomod2.data@data.species == 1)), names(coord)]

#===================================================================
## Geographic Distance Modelling

# geoDist(p, ...)
# p point locations (presence); two column matrix, data.frame or SpatialPoints* object
# ... you must supply a lonlat= argument(logical), unless p is a SpatialPoints* object and has a
#	valid CRS
# ... you can also supply an additional argument 'a' for absence points (currently ignored.); 
#	argument 'a' should be of the same class as argument 'p'

model.sdm <- geoDist(p=occur, lonlat=TRUE)
# save out the model object
EC_Save(model.sdm, EC_OutfileName(filename="model.object",
                                  id_str=species_algo_str, ext="RData"))

# predict for given climate scenario (over constraint region only)
model.proj <- predict(model.sdm, current.climate.scenario@layers[[1]],
                      mask=TRUE, scale=opt.scale)

# remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# save output
EC_SaveModelProjection(model.proj, projection.name, occur.species, species_algo_str)
