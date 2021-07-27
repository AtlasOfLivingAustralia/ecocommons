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
#' @export model_mars


model_mars <- function(x){ # EC.params, EC.env){
  
  # Define the working directory
  # NOTE: This section commented out in source code
  # directories <- set_directories(x$env)
  
  # Set data for modelling
  response_info <- build_response(x$params)  # species data
  predictor_info <- build_predictor(x$params,response_info)  # environmental data
  
  # Set parameters to perform modelling
  model_algorithm <- 'MARS'
  # specific parameters to run MARS algorithm
  model_options_mars <- build_mars_options(x$params)
  # general parameters to run biomod2 package modelling
  model_options_biomod <- build_biomod_options(x$params, response_info,
                                               model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_mars <- c(model_options_biomod$biomod_eval_method,
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
                                    model_options_biomod, model_options_mars)
  
  
  # Save a variable importance plot (VIP) and the model object
  model_save <- save_model_object(model_compute_sdm,
                                  model_algorithm,
                                  model_options_biomod)
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained(predictor_info,model_compute_sdm,
                                      model_options_biomod, response_info,
                                      x$params, 
                                      model_algorithm)
  
  # 2. projection with constraint
  model_projection <- project_model_current(model_compute_sdm, predictor_info,
                                            model_options_biomod,response_info,
                                            x$params, 
                                            model_algorithm)

# RETURN?? not set yet

} # end model_mars()



#________________________________________________________________
# subfunctions to run model_mars()

build_mars_options <- function(a){
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

