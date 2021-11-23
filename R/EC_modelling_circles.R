#' Function to run a Circles within EcoCommons. It is a geographical
#' model to predict a species within a set radius around occurrence records.
#' 
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param dataset_info Dataset object; a nested named list created on EC_build_dataset
#' @param predictor_info Predictor object; a nested named list created on EC_build_predictor
#' @param model_compute Model compute; a bimod2 object
#'
#' @importFrom dismo circles
#'
#' @export EC_modelling_circles

EC_modelling_circles<- function( EC.params, 
                                 response_info,  # from build_response
                                 dataset_info,  # from build_dataset
                                 predictor_info) {   # from EC_build_predictor

  # Set parameters to perform modelling
  model_algorithm <- 'circles'
  # Specific parameters to run Circles algorithm
  opt.d = EC.params$d  # radius around circles in meters; if not specified it is computed from the mean inter-point distance 

  # General parameters to run biomod2 package modelling on Geographical algorithms
  model_options_geographical <- EC_options_geographical (EC.params,
                                                         response_info,
                                                         model_algorithm)
  # Parameters to perform modelling
  if (is.null(opt.d)) {
    model_sdm = dismo::circles(p= dataset_info$occur, lonlat=TRUE)
  } else {
    model_sdm = dismo::circles(p= dataset_info$occur, 
                               d= opt.d, lonlat=TRUE)
  }
  

  #model_save <- EC_save_geographical_model (model_sdm,
  #                                          model_algorithm,
  #                                          model_options_geographical)
  #return(model_sdm)
}
