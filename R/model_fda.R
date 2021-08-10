#' Function to run a Flexible Discriminant Analysis (FDA) within EcoCommons.
#' It is an statistic regression, a classification model based on a mixture
#' of linear regression models, which uses optimal scoring to transform the
#' response variable so that the data are in a better form for linear separation,
#' and multiple adaptive regression splines to generate the discriminant surface
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
#' @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
#' @param EC.env Environment object; as for EC.params. Consider converting to S3
#' 
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' 
#' @export model_fda


model_fda <- function(x){ # EC.params, EC.env){
  
  # Set parameters to perform modelling
  model_algorithm <- 'FDA'
  # specific parameters to run FDA algorithm
  model_options_fda <- build_fda_options(x$params)
  # general parameters to run biomod2 package modelling
  model_options_biomod <- build_biomod_options(x$params, response_info,
                                               model_algorithm )
  
  # Model accuracy statistics
  model_accuracy_fda <- c(model_options_biomod$biomod_eval_method,
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
                                    model_options_biomod, model_options_fda)
  
  
  # Save a variable importance plot (VIP) and the model object
  model_save <- save_model_object(model_compute_sdm,
                                  model_algorithm ,
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
  
} # end model_fda()



#________________________________________________________________
# subfunctions to run model_fda()

build_fda_options<-function(a){
  # Set specific parameters to run FDA algorithm
  list(
    model.options.fda <- list(
      method = EC.params$method  # regression method used in optimal scaling
    )
  )
}