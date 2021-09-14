#' Function to run a Generalized Linear Model (GLM) within EcoCommons. It is an
#' statistic regression, fitted with maximum likelihood estimation for data 
#' with non-normal distribution
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
#' @export EC_modelling_glm


EC_modelling_glm <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset 
  # Set parameters to perform modelling
  model_algorithm <- 'GLM'
  # specific parameters to run GLM algorithm
  model_options_glm <- EC_options_glm (a)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (a, response_info,
                                               model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_glm <- c(model_options_algorithm$biomod_eval_method,
                          model_options_algorithm$dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
    }

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (GLM = model_options_glm)
    
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_algorithm (predictor_info, pa_number_point, a,
                                        model_options_algorithm, model_options_glm)


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
  
} # end EC_modelling_glm()



#________________________________________________________________
# subfunctions to run EC_modelling_glm()

EC_options_glm <- function(a){
  # Set specific parameters to run GLM algorithm
  list(
  	type = a$type,	 # "simple", "quadratic" or "polynomial"; off if myFormula is not NULL
  	interaction_level = a$interaction_level,  # integer; off if myFormula is not NULL
  	# myFormula = NULL,  # specific formula; if not NULL, type and interaction.level are args are switched off
  	test = a$test,  # "AIC", "BIC" or "none"
  	family = a$family,  # e.g. "binomial", "gaussian", "gamma", "inverse.gaussian"
  	mustart = a$mustart,  # starting values for the vector of means
  	control = list(
  		epsilon = a$control_epsilon,  # positive convergence tolerance
  		maxit = a$control_maxit,  # integer giving the maximal number of IWLS iterations
  		trace = a$control_trace  # logical indicating if output should be produced for each iteration
  	)
  )
}

