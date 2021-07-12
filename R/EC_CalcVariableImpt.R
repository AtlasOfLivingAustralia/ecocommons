#' Get the enviromental variables and values used to create the model and
#' generate an output
#'
#' @param out.model 
#' @param model.name 
#' @param num_samples 
#' @param species_algo_str 
#'
#' @return
#' @export EC_CalcVariableImpt
#' 
#' 
EC_CalcVariableImpt <- function(out.model, model.name, num_samples, species_algo_str) {

  if (model.name == "brt") {
    model.values <- matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars <- out.model$var.names
    colnames(model.values) = env.vars
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values <- rbind(out.model@presence, out.model@absence)
    env.vars <- colnames(model.values)
  } else {
    model.values <- out.model@presence
    env.vars <- colnames(model.values)
  }

  if (!(length(model.values)==0)) {

    if (model.name == "brt") {  # predict using actual values
      actual.predictions <- predict(out.model, as.data.frame(model.values),
                                   n.trees = out.model$gbm.call$best.trees, type="response")
    } else {
      actual.predictions = predict(out.model, model.values)
    }
    t
    varimpt.out <- matrix(NA, nrow=length(env.vars), ncol=num_samples+2)  # table to hold the output
    dimnames(varimpt.out) <- list(env.vars, c(paste("sample_", c(1:num_samples, "mean")), "percent"))
    sample.data = model.values  # create a copy of the env data matrix
    
    for (p in 1:ncol(sample.data)) {  # for each predictor variable
      for (s in 1:num_samples) {  # for each num_sample
        
        sample.data[,p] = sample(x=sample.data[,p], replace=FALSE)
        if (model.name == "brt") {
          new.predictions = predict(out.model, as.data.frame(sample.data),
                                    n.trees = out.model$gbm.call$best.trees, type = "response")
        } else {
          new.predictions = predict(out.model, sample.data)
        }
        # calculate correlation between original predictions and new predictions
        varimpt.out[p,s] = 1-max(round(cor(x=actual.predictions,
                                           y=new.predictions, use="pairwise.complete.obs",
                                           method="pearson"), digits=3),0)
      }
    }
    
    # calculate mean variable importance, normalize to percentages, and write results
    varimpt.out[,num_samples+1] = round(rowMeans(varimpt.out, na.rm=TRUE), digits=3)
    varimpt.out[,num_samples+2] = round((varimpt.out[,num_samples+1]/sum(varimpt.out[,num_samples+1]))*100, digits=0)
    EC_WriteCSV(varimpt.out, name=sprintf("biomod2_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate variable importance for ", 
                model.name, "object", sep=" "), stdout())
  }
}
