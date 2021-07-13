########################################################
###   EcoCommons script to BioClim - using biomod2   ###
########################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs BioClim profile, that defines a multidimensional space 
##  bounded by the min and max values of env variables for all occurrences
##  as the species potential range. It uses  the biomod2 package on R to set 
##  up parameters, dismo functions, and the input data generated on the 
##  EcoCommons platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

#===================================================================
## Load and edit your dataset with parameters specific to BioClim

opt.tails <- EC.params$tails # default "both"; use to ignore the left or right tail of the percentile distribution
opt.ext <- NULL #an optional extent object to limit the prediction to a sub-region of 'x'
projection.name <- "current"
species_algo_str <- ifelse(is.null(EC.params$subset), 
                          sprintf("%s_bioclim", occur.species), 
                          sprintf("%s_bioclim_%s", occur.species,
                                  EC.params$subset))

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
  pa_number_point = floor(pa_ratio * nrow(occur))
}

# Format the data as in biomod2. This will also generate the psedo absence points.
biomod2.data <- EC_FormatDataBIOMOD2(true.absen   = absen,
                                  pseudo.absen.points    = pa_number_point,
                                  pseudo.absen.strategy  = EC.params$pa_strategy,
                                  pseudo.absen.disk.min  = EC.params$pa_disk_min,
                                  pseudo.absen.disk.max  = EC.params$pa_disk_max,
                                  pseudo.absen.sre.quant = EC.params$pa_sre_quant,
                                  climate.data           = current.climate.scenario,
                                  occur                  = occur,
                                  species.name           = occur.species,
                                  species_algo_str       = species_algo_str)


#===================================================================
## BioClim Profile


# bioclim(x, p, ...)
# x is a Raster* object or matrix
# p is a two column matrix or SpatialPoints* object
# if p is missing, x is a matrix of values of env vars at known locations of occurrence
# if p is present, it is the location of occurrence and used to extract values for env vars from x,
#       a Raster* object
# NOTE: env vars must be numerical

if (!all(enviro.data.type=="continuous")) {
    stop("bioclim not run because categorical data cannot be used")
} else {
    # run bioclim with matrix of enviro data
    model.sdm <- bioclim(x=occur[,names(current.climate.scenario), drop=FALSE])
    # save out the model object
    EC_Save(model.sdm, EC_OutfileName(filename="model.object",
                                      id_str=species_algo_str, ext="RData"))

    # Do projection over current climate scenario without constraint
    if (enviro.data.genUnconstraintMap && 
       (!is.null(enviro.data.constraints) || enviro.data.generateCHull)) {
        model.proj <- predict(model.sdm, current.climate.scenario.orig, tails=opt.tails)
    
        # remove the current.climate.scenario to release disk space
        EC_RevRasterObject(current.climate.scenario.orig)

        # save output
        EC_SaveModelProjection(model.proj, projection.name, occur.species,
                               species_algo_str, filename_ext="unconstrained")
    }

    # predict for given climate scenario with region constraint
    model.proj <- predict(model.sdm, current.climate.scenario, tails=opt.tails)

    # remove the current.climate.scenario to release disk space
    EC_RevRasterObject(current.climate.scenario)

    # save output
    EC_SaveModelProjection(model.proj, projection.name, occur.species, species_algo_str)

    # evaluate model
    if (!is.null(absen)) {
        EC_SaveDISMOModelEval('bioclim', model.sdm, occur, absen, occur.species)
    }
}
