#' Function to run a Multivariate Adaptive Regression Splines (MARS) within 
#' EcoCommons.It is an statistic regression that builds Builds multiple linear 
#' regression models by partintioning the data and running linear regressions 
#' on each partition
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
#' @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
#' @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' @export EC_modelling_mars


EC_modelling_mars <- function( a,# EC.params
                               response_info,  # from EC_build_response
                               predictor_info,  # from EC_build_predictor
                               dataset_info) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'MARS'
  # specific parameters to run MARS algorithm
  model_options_mars <- EC_options_mars (a)
  # general parameters to run biomod2 package modelling
  model_options_statregr <- EC_options_statregr (a, response_info,
                                                 model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_mars <- c(model_options_statregr$biomod_eval_method,
                          model_options_statregr$dismo_eval_method)
  
  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (MARS = model_options_mars)
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_statregr (predictor_info, pa_number_point, a,
                                        model_options_statregr, model_options_mars)
  
  
  
  # Save a variable importance plot (VIP) and the model object
  model_save <- EC_save_statregr_model (model_compute, model_algorithm,
                                        model_options_statregr)
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained (predictor_info,model_compute,
                                       model_options_statregr, response_info,
                                       a, model_algorithm)
  
  # 2. projection with constraint
  model_projection <- project_model_current (model_compute, predictor_info,
                                             model_options_statregr, response_info,
                                             a, model_algorithm)

# RETURN?? not set yet

} # end EC_modelling_mars()



#________________________________________________________________
# subfunctions to run EC_modelling_mars()

EC_options_mars <- function(a){
  # Set specific parameters to run GAM algorithm
  list(
    type = a$type, #"simple", "quadratic" or "polynomial"; switched off if myFormula is not NULL
    interaction.level = a$interaction_level,   # integer; interaction between variables. Switched off if myFormula is not NULL
    # myFormula = NULL, #specific formula; if not NULL, type and interaction.level args are switched off (N/A for EcoCommons)
    nk =  a$nk, # maximum number of model terms
    penalty= a$penalty, #generalized cross validation (gcv) penalty per knot
    thresh = a$thresh, #forward stepwise stopping threshold
    nprune = a$nprune, #maximum number of terms in the pruned model
    pmethod = a$pmethod #pruning method
  )
}

