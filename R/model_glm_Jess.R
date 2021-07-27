#' Function to run a Generalized Linear Model (GLM) within EcoCommons. It is an
#' statistic regression, fitted with maximum likelihood estimation for data 
#' with non-normal distribution
#' 
#' @param x a list created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' @export model_glm


model_glm <- function(x){ # EC.params, EC.env){

  # Define the working directory
  # NOTE: This section commented out in source code
  # directories <- set_directories(x$env)

  # Set data for modelling
  response_info <- build_response(x$params)  # species data
  predictor_info <- build_predictor(x$params,response_info)  # environmental data

  # Set parameters to perform modelling
  model_algorithm <- 'GLM'
  # specific parameters to run GLM algorithm
  model_options_glm <- build_glm_options(x$params)
  # general parameters to run biomod2 package modelling
  model_options_biomod <- build_biomod_options(x$params, response_info,
                                               model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_glm <- c(model_options_biomod$biomod_eval_method,
                        model_options_biomod$dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- EC.params$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
    }
  
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute_sdm <- mode_compute(predictor_info,
                                   pa_number_point = pa_number_point, x$params,
                                   model_options_biomod, model_options_glm)


  # Save a variable importance plot (VIP) and the model object
  model_save <- save_model_object(model_compute_sdm,
                                model_algorithm,
                                model_options_biomod)
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <- project_model_current_constrained(predictor_info,
                                                                    model_compute_sdm,
                                                                    model_options_biomod,
                                                                    response_info,
                                                                    x$params,
                                                                    model_algorithm)
  
  # 2. projection with constraint
  model_projection <- project_model_current(model_compute_sdm, predictor_info,
                                           model_options_biomod,response_info,
                                           x$params, 
                                           model_algorithm)
  # RETURN?? not set yet
  
} # end model_glm()



#________________________________________________________________
# subfunctions to run model_glm()

build_glm_options <- function(a){
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

