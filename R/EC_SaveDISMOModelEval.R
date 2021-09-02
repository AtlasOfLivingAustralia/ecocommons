#' Save a model evaluation using the dismo R package
#'
#' @param model.name 
#' @param model.obj 
#' @param occur 
#' @param bkgd 
#' @param species.name 
#'
#' @export EC_SaveDISMOModelEval
#' @importFrom dismo evaluate

EC_SaveDISMOModelEval <- function(model.name, model.obj, occur, bkgd, species.name) {
  species_algo_str = paste(species.name, model.name, sep="_")

  # evaluate model using dismo's evaluate
  if (model.name == "brt") {
    model.eval = dismo::evaluate(p=occur, a=bkgd, model=model.obj, n.trees=model.obj$gbm.call$best.trees, type="response")
  } else {
    model.eval = dismo::evaluate(p=occur, a=bkgd, model=model.obj)
  }
  # need predictions and observed values to create confusion matrices for accuracy statistics
  model.fit = c(model.eval@presence, model.eval@absence)
  model.obs = c(rep(1, length(model.eval@presence)), rep(0, length(model.eval@absence)))

  # Call the evaluation script
  res = EC_Performance2D(model.obs, model.fit, species_algo_str, make.plot="dismo", kill.plot=F)
  EC_SaveModelEval(res$performance, res$stats, res$loss.summary, species_algo_str)

  # Create response curves
  EC_CreateResponseCurve(model.obj, model.name, species_algo_str)

  # Calculate variable importance (like biomod2, using correlations between predictions)
  EC_CalcVariableImpt(model.obj, model.name, 3, species_algo_str)

  # Calculate variable importance (like maxent, using decrease in AUC)
  EC_CalcaPermutationVarImpt(model.obj, model.eval, model.name, occur, bkgd, species_algo_str)

  # Create HTML file with accuracy measures
  # bccvl.generateHTML()
}
