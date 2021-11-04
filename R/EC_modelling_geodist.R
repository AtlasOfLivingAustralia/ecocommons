#' Function to run Geographic Distance within EcoCommons. It is a geographical
#' model to predict a species by assuming that a species is more likely to be
#' found closer to an occurrence record. It uses the biomod2 package on R to set
#' up parameters, dismo functions, and the input data generated on the EcoCommons
#' platform.
#'
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param model_compute Model compute; a bimod2 object
#'
#' @importFrom dismo geoDist
#'
#' @export EC_modelling_geodist

EC_modelling_geodist <- function (EC.params,
                                  response_info, # from build_response
                                  model_compute) {

  # Set parameters to perform modelling
  model_algorithm <- 'geodist'
  
  # General parameters to run biomod2 package modelling on Geographical algorithms
  model_options_geographical <- EC_options_geographical (EC.params, response_info,
                                                         model_algorithm)

  # Parameters to perform modelling
  coord <- model_compute@coord
  occur <- coord[c(which(model_compute@data.species == 1)), names(coord)]
  model_sdm <- dismo::geoDist(p= occur, lonlat=TRUE)

  # Save out the model object
  model_save <- EC_save_geographical_model (model_compute,
                                            model_algorithm,
                                            model_options_geographical)
}

