#' Function to run a Convex Hull within EcoCommons. It is a geographical
#' model to predict a species within a minimum bounds of occurrence records.
#' It uses the biomod2 package on R to set up parameters, dismo functions,
#' and the input data generated on the EcoCommons platform.
#'
#' @param a List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param predictor_info Predictor object; a nested named list created on EC_build_predictor
#' @param dataset_info Dataset object; a nested named list created on EC_build_dataset
#' 
#' @importFrom dismo convHull
#'
#' @export EC_modelling_convhull

EC_modelling_convhull<- function( a,# EC.params
                                  response_info,  # from EC_build_response
                                  predictor_info,  # from EC_build_predictor
                                  dataset_info) {  # from EC_build_dataset 

  # Set parameters to perform modelling
  model_algorithm <- 'convhull'
  ## general parameters to run biomod2 package modelling on Geographical algorithms
  model_options_geographical <- EC_options_geographical (a, response_info,
                                                         model_algorithm)

  # Model accuracy statistics
  model_accuracy_convhull <- c(model_options_geographical$biomod_eval_method,
                               model_options_geographical$dismo_eval_method)

  # convhull algorithm does not need pseudo-absence points
  pa_number_point <- 0

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_geographical (predictor_info, a,
                                            model_options_geographical)


  # Parameters to perform modelling
  model_sdm <- dismo::convHull(p=occur)

  # Save out the model object
  model_save <- EC_save_geographical_model (model_compute,
                                            model_algorithm,
                                            model_options_geographical)
}
