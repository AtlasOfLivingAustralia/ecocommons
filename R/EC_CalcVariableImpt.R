EC_CalcVariableImpt <- function(out.model, model.name, num_samples, species_algo_str) {
  # EMG num_samples should be same as biomod.VarImport arg set in
  # 01.init.args.model.current.R

  # get the enviromental variables and values used to create the model
  # EMG this is duplicated from above, should be able to combine
  if (model.name == "brt") {
    model.values = matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars = out.model$var.names
    colnames(model.values) = env.vars
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }

  if (!(length(model.values)==0)) {
    # predict using actual values
    if (model.name == "brt") {
      actual.predictions = predict(out.model, as.data.frame(model.values), n.trees = out.model$gbm.call$best.trees, type="response")
    } else {
      actual.predictions = predict(out.model, model.values)
    }
    # create a table to hold the output
    varimpt.out = matrix(NA, nrow=length(env.vars), ncol=num_samples+2)
    dimnames(varimpt.out) = list(env.vars, c(paste("sample_", c(1:num_samples, "mean")), "percent"))
    # create a copy of the env data matrix
    sample.data = model.values
    # for each predictor variable
    for (p in 1:ncol(sample.data)) {
      # for each num_sample
      for (s in 1:num_samples) {
        # resample from that variables' values, keeping other variable values the same, and predict suitability
        sample.data[,p] = sample(x=sample.data[,p], replace=FALSE)
        # predict using sampled values
        if (model.name == "brt") {
          new.predictions = predict(out.model, as.data.frame(sample.data), n.trees = out.model$gbm.call$best.trees, type = "response")
        } else {
          new.predictions = predict(out.model, sample.data)
        }
        # calculate correlation between original predictions and new predictions
        varimpt.out[p,s] = 1-max(round(cor(x=actual.predictions, y=new.predictions, use="pairwise.complete.obs", method="pearson"), digits=3),0)
      }
    }
    # calculate mean variable importance, normalize to percentages, and write results
    varimpt.out[,num_samples+1] = round(rowMeans(varimpt.out, na.rm=TRUE), digits=3)
    varimpt.out[,num_samples+2] = round((varimpt.out[,num_samples+1]/sum(varimpt.out[,num_samples+1]))*100, digits=0)
    EC_WriteCSV(varimpt.out, name=sprintf("biomod2_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate variable importance for ", model.name, "object", sep=" "), stdout())
  }
}
