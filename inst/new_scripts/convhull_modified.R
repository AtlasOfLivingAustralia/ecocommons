############################################################
###   EcoCommons script to Convex Hull - using biomod2   ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs Convex Hull to predict a species within a minimum 
##  bounds of occurrence records. It is a geographical method for Species 
##  Distribution Modelling. It uses the biomod2 package on R to set up 
##  parameters, dismo functions, and the input data generated on the  
##  EcoCommons platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## Load and edit your dataset with parameters specific to Convex Hull

opt.tails <- EC.params$tails # default "both"; to ignore the left or right tail of the percentile distribution
opt.ext <- NULL #an optional extent object to limit the prediction to a sub-region of 'x'
projection.name <- "current"
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_convhull", occur.species), 
                          sprintf("%s_convhull_%s", occur.species, EC.params$subset))


# model accuracy statistics
# these are available from dismo::evaluate.R NOT originally implemented in biomod2::Evaluate.models.R
dismo.eval.method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")
# and vice versa
biomod.models.eval.meth <- c("KAPPA", "TSS", "ROC", "FAR", "SR", "ACCURACY",
                             "BIAS", "POD", "CSI", "ETS")
# model accuracy statistics - combine stats from dismo and biomod2 for consistent output
model.accuracy <- c(dismo.eval.method, biomod.models.eval.meth)

# Format the data as in biomod2. This will also generate the psedo absence points.
biomod2.data <- bccvl.biomod2.formatData(true.absen       = absen,
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
## Convex Hull Modelling

# convHull(p, ...)
# p point locations (presence), two column matrix, data.frame or SpatialPoints* object
# ... you can supply an argument n (>=1) to get n convex hulls around subset of the points
# ... you can also set n=1:x, to get a set of overlapping polygons consisting of 1 to x parts; i.e.
#   the first polygon has 1 part, the second has 2 parts and x has x parts

# run convhull with occurence data.
model.sdm <- convHull(p=occur)
# save out the model object
EC_Save(model.sdm, EC_OutfileName(filename="model.object",
                                  id_str=species_algo_str, ext="RData"))

# predict for given climate scenario
model.proj <- predict(model.sdm,
                      current.climate.scenario@layers[[1]], mask=TRUE)

# remove the current.climate.scenario to release disk space
EC_RevRasterObject(current.climate.scenario)

# save output
EC_SaveModelProjection(model.proj, projection.name,
                       occur.species, species_algo_str)
