#' Get the environmental variables and values used to create the model and
#' performs permutations
#' 
#' @param out_model 
#' @param model_eval 
#' @param model_name 
#' @param occur 
#' @param bkgd 
#' @param species_algo_str 
#'
#' @export EC_calc_permutation
#' @importFrom dismo evaluate
#'

EC_calc_permutation <- function (out_model, model_eval, model_name, occur, bkgd,
                                 species_algo_str) {
  
  if (model_name == "brt") {
    model_values <- matrix(out_model$data$x, ncol=length(out_model$var.names))
    env_vars <- out_model$var.names
    colnames(model_values) = env_vars
    
  } else if (model_name %in% c("geoIDW", "voronoiHull")) {
    model_values = rbind(out_model@presence, out_model@absence)
    env_vars = colnames(model_values)
    
  } else {
    model_values = out_model@presence
    env_vars = colnames(model_values)
  }

  if (!(length(model_values) == 0)) {
    # get the occurrence and background environmental data used to evaluate the model
    p.swd = occur
    a.swd = bkgd
    # get the AUC from the original model evaluation
    init.auc <- round(model_eval@auc, digits = 3)
    # create a table to hold the output
    permvarimpt.out <- matrix (NA, nrow = length(env_vars), ncol = 4)
    dimnames(permvarimpt.out) <- list(env_vars,
                                      c("init.auc", "sample.auc",
                                        "change.auc", "percent"))
    
    permvarimpt.out[,"init.auc"] = rep(init.auc, length(env_vars))
    # create a copy of the occurrence and background environmental data
    sample.p = p.swd[,env_vars, drop = FALSE]
    sample.a = a.swd[,env_vars, drop = FALSE]
    # check for and remove any NA's present in the data
    no.na.sample.p <- na.omit(sample.p);
    no.na.sample.a <- na.omit(sample.a)
    if (nrow(no.na.sample.p) != nrow(sample.p)) {
      write(paste("EC.calculatePermutationVarImpt(): 
                  NA's were removed from presence data!"), stdout())
    }
    if (nrow(no.na.sample.a) != nrow(sample.a)) {
      write(paste("EC.calculatePermutationVarImpt():
                  NA's were removed from absence data!"), stdout())
    }
    # for each predictor variable
    for (v in 1:length(env_vars)) {
      # resample from that variables' values, keeping other variable values the same
      no.na.sample.p[,v] = sample(x = no.na.sample.p[,v], replace = FALSE)
      no.na.sample.a[,v] = sample(x = no.na.sample.a[,v], replace = FALSE)
      # re-evaluate model with sampled env values
      if (model_name == "brt") {
        sample.eval <- dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a,
                                      model=out_model,
                                      n.trees=out_model$gbm.call$best.trees,
                                      type="response")
      } else {
        sample.eval <- dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a,
                                       model=out_model)
      }
      # get the new auc
      permvarimpt.out[v,"sample.auc"] = round(sample.eval@auc, digits = 3)
    }
    # calculate the difference in auc, normalize to percentages, and write results
    permvarimpt.out[,"change.auc"] = permvarimpt.out[,"init.auc"] - 
      permvarimpt.out[,"sample.auc"]
    
    for (r in 1:nrow(permvarimpt.out)) {
      if (permvarimpt.out[r,"change.auc"] < 0) {  # EMG what if AUC increases?
        permvarimpt.out[r,"change.auc"] = 0
      }
    }
    permvarimpt.out[,"percent"] = round((permvarimpt.out[,"change.auc"] / 
                                           sum(permvarimpt.out[,"change.auc"])) * 100, digits = 0)
    
    EC_write_csv(permvarimpt.out, name = sprintf("maxent_like_VariableImportance_%s.csv", 
                                                 species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate maxent-like variable importance for ",
                model_name, "object", sep=" "), stdout())
  }
}
