#' Get the environmental variables and values used to create the model and
#' performs permutations
#' 
#' @param out.model 
#' @param model.eval 
#' @param model.name 
#' @param occur 
#' @param bkgd 
#' @param species_algo_str 
#'
#' @export EC_CalcPermutationVarImpt
#' @ImportFrom dismo evaluate
#'
EC_CalcPermutationVarImpt <- function(out.model, model.eval,
                                       model.name, occur, bkgd, species_algo_str) {
  if (model.name == "brt") {
    model.values <- matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars <- out.model$var.names
    colnames(model.values) = env.vars
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }

  if (!(length(model.values)==0)) {
    # get the occurrence and background environmental data used to evaluate the model
    p.swd=occur
    a.swd=bkgd
    # get the AUC from the original model evaluation
    init.auc <- round(model.eval@auc, digits=3)
    # create a table to hold the output
    permvarimpt.out <- matrix(NA, nrow=length(env.vars), ncol=4)
    dimnames(permvarimpt.out) <- list(env.vars, c("init.auc", "sample.auc", "change.auc", "percent"))
    permvarimpt.out[,"init.auc"] = rep(init.auc, length(env.vars))
    # create a copy of the occurrence and background environmental data
    sample.p = p.swd[,env.vars, drop=FALSE]
    sample.a = a.swd[,env.vars, drop=FALSE]
    # check for and remove any NA's present in the data
    no.na.sample.p <- na.omit(sample.p);
    no.na.sample.a <- na.omit(sample.a)
    if (nrow(no.na.sample.p) != nrow(sample.p)) {
      write(paste("EC.calculatePermutationVarImpt(): NA's were removed from presence data!"), stdout())
    }
    if (nrow(no.na.sample.a) != nrow(sample.a)) {
      write(paste("EC.calculatePermutationVarImpt(): NA's were removed from absence data!"), stdout())
    }
    # for each predictor variable
    for (v in 1:length(env.vars)) {
      # resample from that variables' values, keeping other variable values the same
      no.na.sample.p[,v] = sample(x=no.na.sample.p[,v], replace=FALSE)
      no.na.sample.a[,v] = sample(x=no.na.sample.a[,v], replace=FALSE)
      # re-evaluate model with sampled env values
      if (model.name == "brt") {
        sample.eval <- dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a,
                                      model=out.model, n.trees=out.model$gbm.call$best.trees,
                                      type="response")
      } else {
        sample.eval <- dismo::evaluate(p=no.na.sample.p, a=no.na.sample.a, model=out.model)
      }
      # get the new auc
      permvarimpt.out[v,"sample.auc"] = round(sample.eval@auc, digits=3)
    }
    # calculate the difference in auc, normalize to percentages, and write results
    permvarimpt.out[,"change.auc"] = permvarimpt.out[,"init.auc"] - permvarimpt.out[,"sample.auc"]
    for (r in 1:nrow(permvarimpt.out)) {
      if (permvarimpt.out[r,"change.auc"] < 0) {  # EMG what if AUC increases?
        permvarimpt.out[r,"change.auc"] = 0
      }
    }
    permvarimpt.out[,"percent"] = round((permvarimpt.out[,"change.auc"]/sum(permvarimpt.out[,"change.auc"]))*100, digits=0)
    EC_WriteCSV(permvarimpt.out, name=sprintf("maxent_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate maxent-like variable importance for ", model.name, "object", sep=" "), stdout())
  }
}
