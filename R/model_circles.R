#' Function to run a Circles within EcoCommons. It is a geographical
#' model to predict a species within a set radius around occurrence records.
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom dismo circles
#'
#' @export model_circles
#' 
#'

model_circles<- function(x){ # EC.params, EC.env){

  # Set parameters to perform modelling
  model_algorithm <- 'circles'
  # Specific parameters to run Circles algorithm
  model_options_circles <- build_circles_options(x$params)
  # General parameters to run biomod2 package modelling on Geographical algorithms
  model_options_biomod <- build_geographical_options(x$params, response_info,
                                               model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_circles <- c(model_options_biomod$biomod_eval_method,
                          model_options_biomod$dismo_eval_method)
  
  # Circles algorithm does not need pseudo-absence points
  pa_number_point <- 0

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute_sdm <- model_geographical_data (predictor_info,
                                    pa_number_point = pa_number_point, x$params,
                                    model_options_biomod)
  

  # Parameters to perform modelling
  if (is.null(opt.d)) {
    model.sdm = dismo::circles(p=occur, lonlat=TRUE)
  } else {
    model.sdm = dismo::circles(p=occur, d=opt.d, lonlat=TRUE)
  }
  
  # Predict for given climate scenario
  model.proj <- predict(model_sdm, predictor_info$current_climate[[1]],
                        mask=TRUE)
  # Save out the model object
  model_save <- save_geographical_model_object(model_compute_sdm,
                                  model_algorithm,
                                  model_options_biomod)
}





#________________________________________________________________
# subfunctions to run model_circles()

build_circles_options <- function(a){
  # Set specific parameters to run "circles"
  list(
    opt.d = a$d # radius around circles; if not specified it is computed from the mean inter-point distance 
 )
}




  