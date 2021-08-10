#' Function to run a Convex Hull within EcoCommons. It is a geographical
#' model to predict a species within a minimum bounds of occurrence records.
#' It uses the biomod2 package on R to set up parameters, dismo functions,
#' and the input data generated on the EcoCommons platform.
#'
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom dismo convHull
#'
#' @export model_convhull
#'

model_convhull<- function(x){ # EC.params, EC.env){

  # Set parameters to perform modelling
  model_algorithm <- 'convhull'
  ## general parameters to run biomod2 package modelling on Geographical algorithms
  model_options_biomod <- build_geographical_options(x$params, response_info,
                                                     model_algorithm)

  # Model accuracy statistics
  model_accuracy_convhull <- c(model_options_biomod$biomod_eval_method,
                              model_options_biomod$dismo_eval_method)

  # convhull algorithm does not need pseudo-absence points
  pa_number_point <- 0

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute_sdm <- model_geographical_data (predictor_info,
                                                pa_number_point = pa_number_point, x$params,
                                                model_options_biomod)


  # Parameters to perform modelling
  model_sdm <- dismo::convHull(p=occur)

  # Save out the model object
  model_save <- save_geographical_model_object(model_compute_sdm,
                                               model_algorithm,
                                               model_options_biomod)
}
