#' Function to run an Artificial Neural Network (ANN) within EcoCommons. It is a
#' machine learning model that uses a hidden layer from linear combinations of 
#' predictor variables to predict species occurrence probabilities
#' 
#' @param a List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param predictor_info Predictor object; a nested named list created on EC_build_predictor
#' @param dataset_info Dataset object; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' 
#' @export EC_modelling_ann


EC_modelling_ann <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'ANN'
  # specific parameters to run ANN algorithm
  model_options_ann <- EC_options_ann (a)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (a, response_info,
                                                 model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_ann <- c(model_options_algorithm$biomod_eval_method,
                          model_options_algorithm$dismo_eval_method)
  
  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (ANN = model_options_ann)
  
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_algorithm (predictor_info, pa_number_point, a,
                                        model_options_algorithm, model_options_ann)
  
  
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
  
} # end EC_modelling_ann()
  

#_____________________________________________________________________________
# subfunctions to run EC_modelling_ann()

EC_options_ann <- function(a){
  # Set specific parameters to run ANN algorithm
  list( NbCV <- a$nbcv, #nb of cross validation to find best size and decay parameters
        size <- a$size, #number of units in the hidden layer
        decay <- a$decay, #parameter for weight decay
        rang <- a$rang, #Initial random weights on [-rang, rang]
        maxit <- a$maxit #maximum number of iterations. Default 100
  )
}
