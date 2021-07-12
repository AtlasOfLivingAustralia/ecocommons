EC_CreateResponseCurve <- function(out.model, model.name, species_algo_str) {
  # Get the enviromental variables and values used to create the model
  if (model.name == "brt") {
    model.values = matrix(out.model$data$x, ncol=length(out.model$var.names))
    env.vars = out.model$var.names
  } else if (model.name %in% c("geoIDW", "voronoiHull")) {
    model.values = rbind(out.model@presence, out.model@absence)
    env.vars = colnames(model.values)
  } else {
    model.values = out.model@presence
    env.vars = colnames(model.values)
  }

  if (!(length(model.values)==0)) {
    # Create a matrix to hold average values for each environmental variable
    mean.values = matrix(data = NA, nrow = 100, ncol = length(env.vars))
    colnames(mean.values) = env.vars
    # For each variable, populate the column with the mean value
    for (i in 1:ncol(mean.values)) {
      mean.values[,i] = rep(mean(model.values[,i], na.rm=TRUE), 100)
    }

    # plot 18 response curves per page
    curvesPerPage = 6*3       # No of rows X No of columns
    for (i in 0:((ncol(mean.values)-1)/curvesPerPage)) {
      png(file=file.path(EC.env$outputdir, sprintf("response_curve_%s_p%d.png", species_algo_str, i)), width=700, height=900)
      par(mfrow = c(6,3)) # No of rows X No of columns

      # Allow each environmental variable to vary, keeping other variable values at average, and predict suitability
      rcurves = list()
      for (j in ((i*curvesPerPage + 1):min((i+1)*curvesPerPage, ncol(mean.values)))) {
        range.values = seq(min(model.values[,j], na.rm=TRUE), max(model.values[,j], na.rm=TRUE), length.out=100)
        temp.data = mean.values
        temp.data[,j] = range.values
        if (model.name == "brt") {
          colnames(temp.data) = env.vars
          new.predictions = predict(out.model, as.data.frame(temp.data), n.trees = out.model$gbm.call$best.trees, type="response")
        } else {
          new.predictions = predict(out.model, temp.data)
        }

        # create separate file for each response curve
        save.name = env.vars[j]
        plot(range.values, new.predictions, ylim=c(0,1), xlab="", ylab="", main=save.name, type="l")
        rug(model.values[,j])

        # save the response curve for later use
        df1 = data.frame(range.values, new.predictions)
        names(df1) <- c(save.name, "")
        rcurves[[save.name]] = df1
      }
      dev.off()

      # Save each response curve
      for (k in 1:length(rcurves))
      {
        ename = env.vars[k + i*curvesPerPage]
        png(file=file.path(EC.env$outputdir, sprintf("%s_response_curve_%s.png", ename, species_algo_str)))
        plot(rcurves[[ename]], ylim=c(0,1), xlab="", ylab="", main=ename, type="l")
        rug(model.values[, k + i*curvesPerPage])
        dev.off()
      }
      rcurves = NULL
    }
  } else {
    write(paste(species_algo_str, ": Cannot create response curves from", model.name, "object", sep=" "), stdout())
  }
}
