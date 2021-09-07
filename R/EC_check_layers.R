#' Function to check if the environmental layers used to project the
#' model are the same as the ones used to create the model object
#' 
#' @param model_obj 
#' @param climatelayers 
#' @param climate_filenames 
#'
#' @export EC_check_layers
#' @importFrom raster dropLayer

EC_check_layers <- function(model_obj, climatelayers, climate_filenames) {
  
  message("Checking environmental layers used for projection")
  if (inherits(model_obj, "DistModel")) {
    model_layers <- colnames(model_obj@presence)
    } else if (inherits(model_obj, "gbm")) {
      model_layers <- summary(model_obj)$var
    
    } else if (inherits(model_obj, "BIOMOD.models.out")) {
      model_layers <- model_obj@expl.var.names
      }
  
  pred_layers <- names(climatelayers)
  
  # check if the env layers were in the original model
  if (sum(!(pred_layers %in% model_layers)) > 0 ){
    
    if (sum(!(model_layers %in% pred_layers)) > 0){
      filenames <- lapply(climate_filenames, function(x) sub("^([^.]*).*",
                                                             "\\1", basename(x)))
      indexes <- match(model_layers, filenames)
      
      for (i in indexes){
        if (!is.na(i)){
          pred_layers[i] = model_layers[i]   # use the corresponding layer name in the model
        }
      }
      names(climatelayers) = pred_layers
    }
    
    message("Dropping environmental layers not used in the original model creation...")
    new.predictors = climatelayers
    for (pl in pred_layers) {
      if (!(pl %in% model_layers)) {
        new.predictors = raster::dropLayer(new.predictors, pl)
      }
    }
    return(new.predictors)
  } else {
    return(climatelayers)
  }
}
