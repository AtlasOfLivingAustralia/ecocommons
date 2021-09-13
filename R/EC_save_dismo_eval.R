#' Save a model evaluation using the dismo R package, creating response curves
#' and calculating variable importance
#'
#' @param model_name 
#' @param model_obj 
#' @param occur 
#' @param bkgd 
#' @param species_name 
#'
#' @importFrom dismo evaluate
#' 
#' #' @export EC_save_dismo_eval

EC_save_dismo_eval <- function (model_name,
                                model_obj,
                                occur,
                                bkgd, 
                                species_name) {
  
  species_algo_str = paste(species_name, model_name, sep="_")

  # evaluate model using dismo's evaluate
  if (model_name == "brt") {
    model_eval = dismo::evaluate (p= occur, a= bkgd, model= model_obj, 
                                  n.trees= model_obj$gbm.call$best.trees, 
                                  type= "response")
  } else {
    model_eval = dismo::evaluate (p= occur, a= bkgd, model= model_obj)
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


#_____________________________________________________________________________
## Subfunctions that are restricted to EC_save_dismo_eval(). listed in the order 
## in which they appear within the function

EC_create_response_curve <- function (out_model,
                                      model_name,
                                      species_algo_str) {
  
  if (model_name == "brt") {
    model_values <- matrix (out_model$data$x, ncol=length(out_model$var.names))
    env_vars = out_model$var.names
    
  } else if (model_name %in% c("geoIDW", "voronoiHull")) {
    model_values <- rbind(out_model@presence, out_model@absence)
    env_vars = colnames(model_values)
    
  } else {
    model_values <- out_model@presence
    env_vars <- colnames(model_values)
  }
  
  if (!(length(model_values) == 0)) {
    # Create a matrix to hold average values for each environmental variable
    mean_values <- matrix(data = NA, nrow = 100, ncol = length(env_vars))
    colnames(mean_values) <- env_vars
    
    # For each variable, populate the column with the mean value
    for (i in 1:ncol(mean_values)) {
      mean_values[,i] = rep(mean(model_values[,i], na.rm=TRUE), 100)
    }
    
    # plot 18 response curves per page
    curvesPerPage = 6 * 3  # no of rows X No of columns
    for (i in 0:((ncol(mean_values) - 1) / curvesPerPage)) {
      png(file = file.path (EC.env$outputdir,
                            sprintf("response_curve_%s_p%d.png",
                                    species_algo_str, i)), width = 700, height = 900)
      
      par(mfrow = c(6,3)) # No of rows X No of columns
      
      # Allow each environmental variable to vary, keeping other variable 
      # at average, and predict suitability
      rcurves = list()
      
      for (j in ((i * curvesPerPage + 1):min((i + 1) * curvesPerPage,
                                             ncol(mean_values)))) {
        range_values = seq(min(model_values[,j], na.rm=TRUE),
                           max(model_values[,j], na.rm=TRUE), length.out=100)
        
        temp_data = mean_values
        
        temp_data[,j] = range_values
        if (model_name == "brt") {
          colnames(temp_data) = env_vars
          new_predictions = predict (out_model, as.data.frame(temp_data),
                                     n.trees = out_model$gbm.call$best.trees, 
                                     type="response")
        } else {
          new_predictions = predict (out_model, temp_data)
        }
        
        # create separate file for each response curve
        save_name = env_vars[j]
        plot(range_values, new_predictions, ylim=c(0,1),
             xlab = "", ylab = "", main = save_name, type="l")
        rug(model_values[,j])
        
        df1 <- data.frame(range_values, new_predictions)  # save the response curve for later use
        names(df1) <- c(save_name, "")
        rcurves[[save_name]] = df1
      }
      dev.off()
      
      for (k in 1:length(rcurves)) # Save each response curve
      {
        ename = env_vars[k + i * curvesPerPage]
        png(file=file.path(EC.env$outputdir,
                           sprintf("%s_response_curve_%s.png", ename, species_algo_str)))
        plot(rcurves[[ename]], ylim = c(0,1), xlab= "", ylab= "", main = ename, type="l")
        rug(model_values[, k + i * curvesPerPage])
        dev.off()
      }
      rcurves = NULL
    }
  } else {
    write(paste(species_algo_str, ": Cannot create response curves from",
                model_name, "object", sep = " "), stdout())
  }
}


#_____________________________________________________________________________


EC_calc_variable_impt <- function(out.model,
                                  model.name,
                                  num_samples, 
                                  species_algo_str) {
  
  if (model.name == "brt") {
    model.values <- matrix(out.model$data$x,
                           ncol=length(out.model$var.names))
    
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
                                    n.trees = out.model$gbm.call$best.trees,
                                    type="response")
    } else {
      actual.predictions = predict(out.model, model.values)
    }
    t
    varimpt.out <- matrix(NA, nrow=length(env.vars), ncol=num_samples+2)  # table to hold the output
    dimnames(varimpt.out) <- list(env.vars, c(paste("sample_",
                                                    c(1:num_samples, "mean")),
                                              "percent"))
    sample.data = model.values  # create a copy of the env data matrix
    
    for (p in 1:ncol(sample.data)) {  # for each predictor variable
      for (s in 1:num_samples) {  # for each num_sample
        
        sample.data[,p] = sample(x=sample.data[,p], replace=FALSE)
        if (model.name == "brt") {
          new.predictions = predict(out.model, as.data.frame(sample.data),
                                    n.trees = out.model$gbm.call$best.trees,
                                    type = "response")
        } else {
          new.predictions = predict(out.model, sample.data)
        }
        # calculate correlation between original predictions and new predictions
        varimpt.out[p,s] = 1-max(round(cor(x=actual.predictions,
                                           y=new.predictions,
                                           use="pairwise.complete.obs",
                                           method="pearson"), digits=3),0)
      }
    }
    
    # calculate mean variable importance, normalize to percentages, and write results
    varimpt.out[,num_samples+1] = round(rowMeans(varimpt.out, na.rm=TRUE), digits=3)
    varimpt.out[,num_samples+2] = round((varimpt.out[,num_samples+1]/sum(varimpt.out[,num_samples+1]))*100, digits=0)
    EC_write_csv(varimpt.out, name=sprintf("biomod2_like_VariableImportance_%s.csv", species_algo_str))
  } else {
    write(paste(species, ": Cannot calculate variable importance for ", 
                model.name, "object", sep=" "), stdout())
  }
}


#_____________________________________________________________________________


EC_calc_permutation <- function (out_model,
                                 model_eval,
                                 model_name,
                                 occur,
                                 bkgd,
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
