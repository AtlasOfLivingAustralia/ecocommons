#' Function to run a Circles within EcoCommons. It is a geographical
#' model to predict a species within a set radius around occurrence records.
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom dismo circles
#'
#' @export EC_modelling_circles
#' 
#'

EC_modelling_circles<- function( a,  # EC.params
                                 response_info,  # from EC_build_response
                                 predictor_info,  # from EC_build_predictor
                                 dataset_info) {  # from EC_build_dataset 

  # Set parameters to perform modelling
  model_algorithm <- 'circles'
  # Specific parameters to run Circles algorithm
  model_options_circles <- EC_options_circles (a)
  # General parameters to run biomod2 package modelling on Geographical algorithms
  model_options_geographical <- EC_options_geographical (a, response_info,
                                                         model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_circles <- c(model_options_geographical$biomod_eval_method,
                              model_options_geographical$dismo_eval_method)

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_geographical (predictor_info, a,
                                            model_options_geographical)
  

  # Parameters to perform modelling
  if (is.null(model_options_circles$opt.d)) {
    model.sdm = dismo::circles(p=dataset_info$occur, lonlat=TRUE)
  } else {
    model.sdm = dismo::circles(p=dataset_info$occur, 
                               d=model_options_circles$opt.d, lonlat=TRUE)
  }
  
  # Predict for given climate scenario
  model.proj <- predict(model_sdm, predictor_info$current_climate[[1]],
                        mask=TRUE)
  # Save out the model object
  model_save <- EC_save_geographical_model (model_compute,
                                            model_algorithm,
                                            model_options_geographical)
}





#________________________________________________________________
# subfunctions to run model_circles()

EC_options_circles <- function(a){
  # Set specific parameters to run "circles"
  list(
    opt.d = a$d # radius around circles; if not specified it is computed from the mean inter-point distance 
 )
}
  