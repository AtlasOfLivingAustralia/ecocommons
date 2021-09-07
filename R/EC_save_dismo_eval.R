#' Save a model evaluation using the dismo R package
#'
#' @param model_name 
#' @param model_obj 
#' @param occur 
#' @param bkgd 
#' @param species_name 
#'
#' @export EC_save_dismo_eval
#' @importFrom dismo evaluate

EC_save_dismo_eval <- function (model_name, model_obj, occur, bkgd, 
                                species_name) {
  
  species_algo_str = paste(species_name, model_name, sep="_")

  # evaluate model using dismo's evaluate
  if (model_name == "brt") {
    model_eval = dismo::evaluate (p=occur, a=bkgd, model=model_obj, 
                                  n.trees=model_obj$gbm.call$best.trees, 
                                  type="response")
  } else {
    model_eval = dismo::evaluate (p=occur, a=bkgd, model=model_obj)
  }
  # predictions and observed values to create confusion matrices for accuracy statistics
  model_fit = c(model_eval@presence, model_eval@absence)
  model_obs = c(rep(1, length(model_eval@presence)), 
                rep(0, length(model_eval@absence)))

  # Call the evaluation script
  res = EC_performance_2D (model_obs, model_fit, species_algo_str, 
                          make.plot="dismo", kill.plot=F)
  
  EC_save_eval (res$performance, res$stats, 
                   res$loss.summary, species_algo_str)

  # Create response curves
  EC_create_response_curve (model_obj, model_name, species_algo_str)

  # Calculate variable importance (like biomod2, using correlations between predictions)
  EC_calc_variable_impt (model_obj, model_name, 3, species_algo_str)

  # Calculate variable importance (like maxent, using decrease in AUC)
  EC_calc_permutation (model_obj, model_eval, model_name, occur,
                       bkgd, species_algo_str)

  # Create HTML file with accuracy measures
  # bccvl.generateHTML()
}
