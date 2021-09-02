#' Title
#'
#' @param loaded.names 
#' @param biomod.model 
#' @param species_algo_str 
#'
#' @importFrom biomod2 get_evaluations
#' @importFrom biomod2 get_predictions
#' @importFrom biomod2 get_formal_data
#' @importFrom biomod2 get_variables_importance
#' @importFrom biomod2 get_formal_data
#' 

EC_SaveBIOMODModelEval <- function(loaded.names, biomod.model, species_algo_str) {

  evaluation = get_evaluations(biomod.model)

  # Get the model predictions and observed values. Predictions is a 4-dimensional array (Predictions, Algorithm, Model run, PseudoAbsence Run)
  predictions = get_predictions(biomod.model)
  total_models = length(dimnames(predictions)[[3]])

  obs = get_formal_data(biomod.model, "resp.var")
  # In case of pseudo-absences we might have NA values in obs so replace them with 0
  obs = replace(obs, is.na(obs), 0)

  for ( i in 1:total_models )
  {
    model_name = dimnames(predictions)[[3]][i]  # will be FULL or RUN1 for eg
    model_predictions = predictions[,,i,]

    if (sum(is.na(model_predictions)) == length(model_predictions))
    {
      # Warn that model n is being ignored. It most probably failed to build.
      warning(sprintf("Warning: Model %i failed to generate. Not generating stats", i), immediate.=T)
      next
    }

    res = EC_Performance2D(obs, model_predictions / 1000, species_algo_str, make.plot=model_name, kill.plot=F)
    EC_SaveModelEval(res$performance, res$stats, res$loss.summary, species_algo_str)

    # get and save the variable importance estimates
    variableImpt = get_variables_importance(biomod.model)
    if (!is.na(variableImpt)) {
      #EMG Note this will throw a warning message if variables (array) are returned
      EC_write_csv(variableImpt,
                  name=paste("variableImportance", model_name, species_algo_str, "csv", sep="."))
    } else {
      message("VarImport argument not specified during model creation!")
      #EMG must create the model with the arg "VarImport" != 0
    }
  }

  # save response curves (Elith et al 2005)
  for(name in loaded.names)
  {
    env_data = get_formal_data(biomod.model,"expl.var")
    png(file=file.path(EC.env$outputdir, sprintf("mean_response_curves_%s.png", name)))
    test <- response.plot2(models = name,
                           Data = env_data,
                           show.variables = get_formal_data(biomod.model,"expl.var.names"),
                           fixed.var.metric = "mean")
    dev.off()

    # save individual response curves
    for (envname in names(test))
    {
      png(file=file.path(EC.env$outputdir, sprintf("%s_mean_response_curves_%s.png", envname, name)))
      plot(test[[envname]], type='l', ylim=c(0, 1.0), main=envname, xlab="", ylab="")
      rug(env_data[[envname]])
      dev.off()
    }
  }
}
