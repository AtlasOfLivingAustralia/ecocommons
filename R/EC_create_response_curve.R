#' Create a Response Curve based on the environmental variables and 
#' and model results
#'
#' @param out_model 
#' @param model_name 
#' @param species_algo_str 
#'
#' @export EC_create_response_curve

EC_create_response_curve <- function (out_model, model_name, species_algo_str) {
  
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
