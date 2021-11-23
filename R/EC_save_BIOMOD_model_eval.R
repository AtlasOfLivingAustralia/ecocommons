#' Save evaluation models for all statistics in csv format
#'
#' @param loaded_model from project_model_current, extract the root of filenames used by biomod2 to save model results
#' @param model_sdm    the final model as result of biomod2
#' @param species_algo_str algorithm being calculated plus the name of the species
#' 
#' @importFrom biomod2 get_evaluations
#' @importFrom biomod2 get_formal_data
#' @importFrom biomod2 get_predictions
#' @importFrom biomod2 get_variables_importance
#' @importFrom biomod2 response.plot2
#'
#' @export EC_save_BIOMOD_model_eval



EC_save_BIOMOD_model_eval <- function(loaded_model,         # as defined by bimod2
                                      model_sdm,            # from model created using biomod2::BIOMOD_Modeling
                                      species_algo_str) {   # from EC_options_algorithm
  
  evaluation = biomod2::get_evaluations(model_sdm)
  
  # Get the model predictions and observed values. 
  predictions  = biomod2::get_predictions(model_sdm)  # Predictions is a 4-dimensional array (Predictions, Algorithm, Model run, PseudoAbsence Run)
  total_models = length(dimnames(predictions)[[3]])
  
  obs = biomod2::get_formal_data(model_sdm, "resp.var")
  obs = replace(obs, is.na(obs), 0) # In case of pseudo-absences we might have NA values in obs so replace them with 0
  
  for ( i in 1:total_models ) {
    model_name = dimnames(predictions)[[3]][i]  # will be FULL or RUN1 for eg
    model_predictions = predictions[,,i,]
    
    if (sum(is.na(model_predictions)) == length(model_predictions)) {
      # Warn that model n is being ignored. It most probably failed to build.
      warning(sprintf("Warning: Model %i failed to generate. Not generating stats", i), immediate.=T)
      next
    }
    
    res = EC_performance_2D (obs              = obs,
                             pred             = model_predictions / 1000,
                             species_algo_str = species_algo_str,
                             make.plot        = model_name,
                             kill.plot        = FALSE)
    
    EC_save_model_eval (res$performance,
                        res$stats,
                        res$loss.summary,
                        res$species_algo_str)
    
    # get and save the variable importance estimates
    variableImpt = biomod2::get_variables_importance (model_sdm)
    if (!is.na (variableImpt)) {
      #EMG Note this will throw a warning message if variables (array) are returned
      EC_write_csv( variableImpt,
                    name = paste("variableImportance", model_name, species_algo_str,
                                "csv", sep = "."))
    } else {
      message ("VarImport argument not specified during model creation!")
    }
  }
  
  # save response curves (Elith et al 2005)
  for(name in loaded_model) {
    
    env_data  = biomod2::get_formal_data(model_sdm, "expl.var")
    png (file = file.path(EC.env$outputdir, sprintf("mean_response_curves_%s.png", name)))
    test <- biomod2::response.plot2 (models           = name,
                                     Data             = env_data,
                                     show.variables   = biomod2::get_formal_data (model_sdm, "expl.var.names"),
                                     fixed.var.metric = "mean")
    dev.off()
    
    # save individual response curves
    for (envname in names(test)) {
      png (file = file.path(EC.env$outputdir, sprintf("%s_mean_response_curves_%s.png",
                                                     envname, name)))
      plot(test[[envname]], type = 'l', ylim = c(0, 1.0), main = envname,
           xlab = "", ylab = "")
      rug(env_data[[envname]])
      dev.off()
    }
  }
}


#_____________________________________________________________________________
# subfunctions to run EC_save_BIOMOD_model_eval()


EC_save_model_eval <- function(out_performance,
                               out_stats,
                               out_lossfunction,
                               species_algo_str){
  
  EC_write_csv(data.frame(out_performance), 
                  name = paste0(paste("Evaluation-data", 
                                      species_algo_str, sep="_"), ".csv"))
  EC_write_csv(data.frame(out_stats), 
                  name = paste0(paste("Evaluation-statistics", 
                                      species_algo_str, sep="_"), ".csv"))
  EC_write_csv(data.frame(out_lossfunction), 
                  name = paste0(paste("Loss-function-intervals-table", 
                                      species_algo_str, sep="_"), ".csv"))
}