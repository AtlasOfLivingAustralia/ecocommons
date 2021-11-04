#' Functions to set up the data to run the algorithms scripts
#'
#' The EcoCommons platform computes many Species Distribution Modelling
#' algorithms to be selected by the user. However, some parameters are repeated
#' across algorithms. This functions aims to compute simple lists for parameters
#' used in all models. It includes response, predictors and constraints,
#' necessary to build the dataset to be used. 
#'
#' @param a List created from a json file, containing source_file$params
#' @param r from response_info Response object; a nested named list created on EC_build_response
#' @param p from predictor_info Predictor object; a nested named list created on EC_build_predictor
#' @param c from constraint_info Constraint object; a nested named list created on EC_build_constraint
#'
#' @keywords internal

EC_build_response <- function(a){  # formerly EC.params
  # Set list with species data information for modelling
  list(
    occur_data    = a$species_occurrence_dataset$filename,
    occur_species = a$species_occurrence_dataset$species,
    month_filter  = a$species_filter,
    absen_data    = a$species_absence_dataset$filename,  # background/pseudo-absence points- 2 column matrix
    pa_ratio      = a$pa_ratio, # pseudo-absence ratio; used later on to create pseudo-absence points
    
    # OPTIONAL bias file from user
    bias_data         = a$bias_dataset$filename
  )
}


#_____________________________________________________________________________


EC_build_predictor <- function(a) { # formerly EC.params
  # Set list with environmental data for modelling
  list(
    current  =  lapply(a$environmental_datasets, function(x) x$filename),
    type     =  lapply(a$environmental_datasets, function(x) x$type),  # continuous or categorical
    layer    =  lapply(a$environmental_datasets, function(x) x$layer)  # current layer names
  )
}


#_____________________________________________________________________________


EC_build_constraint <- function(a) {  # formerly EC.params
  # Read (if exist) a region constraint modelling to run the algorithms, 
  # either specified by the user or applying a Convex Hull; can be NULL
  list (
    seed = a$random_seed,
    #if (is.null(seed)) {
    #  seed = runif(1, -2^31, 2^31-1)
    #}
    #seed = as.integer(seed)
    #set.seed(seed)
    constraints = readLines(a$modelling_region$filename),  # geographic constraints
    #if (!is.null (a$modelling_region$filename)) {
    #  constraints <- readLines(a$modelling_region$filename)
    #  },
    # whether to create convex-hull polygon of occurrence dataset to constraint
    generateCHull = ifelse(is.null(a$generate_convexhull), FALSE,
                            as.logical(a$generate_convexhull)),
    # whether to generate an unconstrained map or not. TRUE by default
    genUnconstraintMap = ifelse(is.null(a$unconstraint_map), TRUE,
                                 as.logical(a$unconstraint_map)),
    # resampling (up / down scaling); if scale_down is TRUE, return 'lowest'
    resampling = ifelse(is.null(a$scale_down) ||
                        as.logical(a$scale_down),
                        'highest', 'lowest')
    )
}


#_____________________________________________________________________________


EC_build_dataset <- function( p,  # formerly predictor_info,
                              c,  # formerly constraint_info 
                              r) {  # formerly response_info  
# read current climate data
current_climate <- EC_raster_stack(p$current, p$type, p$layer,
                                  resamplingflag = c$resampling)
  
# read in the occurrence and absence data (if available); optional month filtering
occur <- EC_read_sp(r$occur_data,
                   r$month_filter)
absen <- EC_read_sp(r$absen_data,
                   r$month_filter)

# if no constrain regions is given, geographically constrained modelling with Convex Hull
constrained_results <- EC_SDM_geoconstrained(current_climate, occur, absen,
                                              c$constraints, c$generateCHull)

if (r$pa_ratio > 0) {
  pa_number_point = floor(r$pa_ratio * nrow(constrained_results$occur))
} else {
  pa_number_point = 0
}

    
  # save a copy of the climate dataset
mylist<- list (
  "current_climate_orig" = current_climate,
  "occur" = constrained_results$occur,
  "absen" = constrained_results$absen,
  "pa_number_point" = pa_number_point,
  "current_climate" = constrained_results$masked_raster
  )

return(mylist)

}
