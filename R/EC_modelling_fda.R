#' Function to run a Flexible Discriminant Analysis (FDA) within EcoCommons.
#' It is an statistic regression, a classification model based on a mixture
#' of linear regression models, which uses optimal scoring to transform the
#' response variable so that the data are in a better form for linear separation,
#' and multiple adaptive regression splines to generate the discriminant surface
#' 
#' @param a a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#' 
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' 
#' @export modelling_fda


modelling_fda <- function( a,# EC.params
                           response_info,  # from EC_build_response
                           predictor_info,  # from EC_build_predictor
                           dataset_info) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'FDA'
  # specific parameters to run FDA algorithm
  model_options_fda <- EC_options_fda (a)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (a, response_info,
                                                 model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_fda <- c(model_options_algorithm$biomod_eval_method,
                          model_options_algorithm$dismo_eval_method)
  
  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (FDA = model_options_fda)
  
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_algorithm (predictor_info, pa_number_point, a,
                                        model_options_algorithm, model_options_fda)
  
  
  # Save a variable importance plot (VIP) and the model object
  model_save <- EC_save_algorithm_model (model_compute, model_algorithm,
                                        model_options_algorithm)
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained (predictor_info,model_compute,
                                       model_options_algorithm, response_info,
                                       a, model_algorithm)
  
  # 2. projection with constraint
  model_projection <- project_model_current (model_compute, predictor_info,
                                            model_options_algorithm, response_info,
                                            a, model_algorithm)

  # RETURN?? not set yet
  
} # end modelling_fda()



#________________________________________________________________
# subfunctions to run modelling_fda()

EC_options_fda <- function(a){
  # Set specific parameters to run FDA algorithm
  list(
    model.options.fda <- list(
      method = a$method  # regression method used in optimal scaling
    )
  )
}