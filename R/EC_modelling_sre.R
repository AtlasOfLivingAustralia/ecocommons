#' Function to run a Surface Range Envelope (SRE) within EcoCommons. It is a
#' profile model to define a multidimensional space bounded by the min and max
#' values of environment variables for all occurrences as the species potential
#' range.
#'
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @export EC_modelling_sre

EC_modelling_sre<- function(x){ # EC.params, EC.env){

# Set parameters to perform modelling
  model_algorithm <- 'SRE'
  # No specific parameter is required to run SRE
  # General parameters to run biomod2 package modelling on Geographical algorithms
  model_options_statregr <- EC_options_statregr (a, response_info,
                                                 model_algorithm)

  # Model accuracy statistics
  model_accuracy_sre <- c(model_options_statregr$biomod_eval_method,
                          model_options_statregr$dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (SRE = model_options_statregr)

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_statregr (predictor_info, pa_number_point, a,
                                        model_options_statregr)
}

