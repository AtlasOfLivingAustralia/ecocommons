#' Function to run Voronoi Hull Model within EcoCommons. It is a geographical
#' model to predict species occur within Voronoi Hulls which surround known
#' occurrences. It uses the biomod2 package on R to set up parameters, dismo
#' functions, and the input data generated on the EcoCommons platform.
#'
#'@param a a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom dismo voronoiHull
#'
#' @export EC_modelling_voronoihull

EC_modelling_voronoihull<- function(  a,# EC.params
                                      response_info,  # from EC_build_response
                                      predictor_info,  # from EC_build_predictor
                                      dataset_info) {  # from EC_build_dataset 

  # Set parameters to perform modelling
  model_algorithm <- 'voronoihull'
  ## general parameters to run biomod2 package modelling on Geographical algorithms
  model_options_geographical <- EC_options_geographical (a, response_info,
                                                         model_algorithm)

  # Model accuracy statistics
  model_accuracy_voronoihull <- c(model_options_geographical$biomod_eval_method,
                                  model_options_geographical$dismo_eval_method)


  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_geographical (predictor_info, a,
                                            model_options_geographical)


  # Parameters to perform modelling
  coord <- model_compute@coord
  occur <- coord[c(which(model_compute@data.species == 1)), names(coord)]
  absen <- coord[c(which(model_compute@data.species == 0 |
                           is.na(model_compute@data.species))), names(coord)]

  model_sdm <- voronoiHull (p = occur, a = absen)

  # Save out the model object
  model_save <- EC_save_geographical_model (model_compute,
                                            model_algorithm,
                                            model_options_geographical)
}

