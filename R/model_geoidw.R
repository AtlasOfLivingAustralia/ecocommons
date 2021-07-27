#' Function to run Inverse Distance Weighted Model within EcoCommons. It is a
#' geographical model to predict a species probabilities in unknown locations
#' by averaging values of nearby occurrences. It uses the biomod2 package on R
#' to set up parameters, dismo functions, and the input data generated on the
#' EcoCommons platform.
#'
#'@param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom dismo geoIDW
#'
#' @export model_geoidw
#'
#'

model_geoidw<- function(x){ # EC.params, EC.env){

  # Define the working directory
  # NOTE: This section commented out in source code
  # directories <- set_directories(x$env)

  # Set data for modelling
  response_info <- build_response(x$params)  # species data
  predictor_info <- build_predictor(x$params,response_info)  # environmental data


  # Set parameters to perform modelling
  model_algorithm <- 'geoidw'
  ## general parameters to run biomod2 package modelling on Geographical algorithms
  model_options_biomod <- build_geographical_options(x$params, response_info,
                                                     model_algorithm)

  # Model accuracy statistics
  model_accuracy_geoidw <- c(model_options_biomod$biomod_eval_method,
                              model_options_biomod$dismo_eval_method)

  # geoidw algorithm does not need pseudo-absence points
  pa_number_point <- 0
  pa_number_point = 0
  if (pa_ratio > 0) {
    pa_number_point = floor(pa_ratio * nrow(occur))
  }

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute_sdm <- model_geographical_data (predictor_info,
                                                pa_number_point = pa_number_point, x$params,
                                                model_options_biomod)


  # Parameters to perform modelling
  coord <- model_compute_sdm@coord
  occur <- coord[c(which(model_compute_sdm@data.species == 1)), names(coord)]
  absen <- coord[c(which(model_compute_sdm@data.species == 0 | is.na(model_compute_sdm@data.species))), names(coord)]

  model_sdm <- geoIDW (p = occur, a = absen)

  # Save out the model object
  model_save <- save_geographical_model_object(model_compute_sdm,
                                               model_algorithm,
                                               model_options_biomod)
}

