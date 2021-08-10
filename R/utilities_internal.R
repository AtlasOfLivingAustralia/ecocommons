#' Define functions to run algorithms scripts
#'
#' The EcoCommons plattform computes many Species Distribution Modelling
#' algorithms to be selected by the user. Each algorithm requires a specific
#' set of parameters. This funtion simply aims to compute simple names for
#' parameters and assign the necessary parameters to run specific algorithm.
#'
#' @param param 
#' @param value 
#'
#' 

# NOTE: set_directories() currently not used
set_directories <- function(a){  # formerly EC.params
  # Define the working directory
  list(
    script = normalizePath(a$scriptdir),
    input =  normalizePath(a$inputdir),
    output = normalizePath(a$outputdir)
  )
}



build_response <- function(a){
  # Set list with species data info for modelling
  list(
    occur_data    = a$species_occurrence_dataset$filename,
    occur_species = a$species_occurrence_dataset$species,
    month_filter  = a$species_filter,
    absen_data    = a$species_absence_dataset$filename,  # background/pseudo-absence points- 2 column matrix
    # OPTIONAL bias file from user
    bias_data     = a$bias_dataset$filename,
    if (is.null (a$bias_dataset$filename)) {
      a$bias_dataset$filename = NULL
  } else {
      bias_data = readLines(a$bias_dataset$filename)
  }
  )
}


    



build_predictor <- function(
  a,  # formerly EC.params
  response_info) { # from build_response()
  # Set environmental data for modelling; first, a list with elements info
  # Then, reading environment data and applying modelling constraint, if necessary
  list(
    current = lapply(a$environmental_datasets, function(x) x$filename),
    type    = lapply(a$environmental_datasets, function(x) x$type),  # continuous or categorical
    layer   = lapply(a$environmental_datasets, function(x) x$layer),  # current layer names
    constraints <- a$modelling_region,  # geographic constraints
    if (is.null (a$modelling_region$filename)) {  # I added this; just if makes sense
       a$modelling_region$filename = NULL
       } else {
         bias_data = readLines(a$modelling_region$filename)
         },
    
    # whether to create convex-hull polygon of occurrence dataset to constraint
    generateCHull <- ifelse(is.null(a$generate_convexhull), FALSE,
                            as.logical(a$generate_convexhull)),
    
    # whether to generate an unconstrained map or not. TRUE by default
    genUnconstraintMap <- ifelse(is.null(a$unconstraint_map), TRUE,
                                 as.logical(a$unconstraint_map)),
    
    # resampling (up / down scaling); if scale_down is TRUE, return 'lowest'
    resampling <- ifelse(is.null(a$scale_down) ||
                           as.logical(a$scale_down),
                         'highest', 'lowest')
  )
  
  # read current climate data
  current_climate <- EC_EnviroStack(current, type, layer,
                                    resamplingflag = resampling)
  
  # read in the necessary observation, background and environmental data
  occur <- EC_ReadSp(response_info$occur_data,
                     response_info$month_filter)
  absen <- EC_ReadSp(response_info$absen_data,
                     response_info$month_filter)
  
  # geographically constrained modelling
  if (!is.null(constraints) || generateCHull) {
    constrained_results <- EC_SDMGeoConstrained(current_climate,occur, absen,
                                                constraints, generateCHull);
    
    # save a copy of the climate dataset
    current_climate_orig <- current_climate  
    current_climate <- constrained_results$raster
    occur <- constrained_results$occur
    absen <- constrained_results$absen
  }
}
