# Function to check if the environmental layers used to project the
# model are the same as the ones used to create the model object

EC_CheckLayers <- function(model.obj, climatelayers, climate_filenames) {
  message("Checking environmental layers used for projection")
  if (inherits(model.obj, "DistModel")) {
    model.layers = colnames(model.obj@presence)  # dismo package
  } else if (inherits(model.obj, "gbm")) {
    model.layers = summary(model.obj)$var      # brt package
  } else if (inherits(model.obj, "BIOMOD.models.out")) {      # biomod package
    model.layers = model.obj@expl.var.names
  }
  
  pred.layers = names(climatelayers)
  
  # check if the env layers were in the original model
  if (sum(!(pred.layers %in% model.layers)) > 0 ){

    if (sum(!(model.layers %in% pred.layers)) > 0){
      filenames = lapply(climate_filenames, function(x) sub("^([^.]*).*", "\\1", basename(x)))
      indexes = match(model.layers, filenames)
      for (i in indexes){
        if (!is.na(i)){
          pred.layers[i] = model.layers[i]   # use the corresponding layer name in the model
        }
      }
      names(climatelayers) = pred.layers
    }
    
    message("Dropping environmental layers not used in the original model creation...")
    new.predictors = climatelayers
    for (pl in pred.layers) {
      if (!(pl %in% model.layers)) {
        new.predictors = dropLayer(new.predictors, pl)
      }
    }
    return(new.predictors)
  } else {
    return(climatelayers)
  }
}
