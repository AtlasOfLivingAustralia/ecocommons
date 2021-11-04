#' Subfunctions for Species Distribution Modelling (SDMs) Geographical 
#' algorithms. It builds evaluation parameters, set data on biomod2
#' structure and save the result
#' 
#' @keywords internal

EC_options_geographical <- function( a, # formerly EC.params
                                     response_info, # from build_response()
                                     model_algorithm){
  # Build evaluation parameters for biomod2 and dismo modelling
  
  list(
    # Model accuracy statistics
    # some are available from biomod2::Evaluate.models.R - "KAPPA", "TSS", "ROC" ,
    #"FAR", "SR", "ACCURACY", "BIAS", "POD", "CSI", "ETS"
    # others from dismo::evaluate.R - "ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR"
    model_accuracy = c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY", "BIAS",
                       "POD", "CSI", "ETS", "ODP", "TNR", "FPR", "FNR", "NPP",
                       "MCR", "OR"),

    opt.scale = ifelse(is.null(a$scale), 1000, a$scale),
    opt.tails = a$tails, # default "both"; use to ignore the left or right
    # tail of the percentile distribution
    opt.ext = NULL, #an optional extent object to limit the prediction to a sub-region of 'x'
    projection_name = "current",
    species_algo_str = if (!is.null(a$subset)) {
                              paste0("_", model_algorithm, "_",
                                     response_info$occur_species)
      }
  )
}


#_____________________________________________________________________________


EC_save_geographical_model<-function (model_compute,  # from model_compute
                                      model_algorithm,
                                      model_options_geographical){  # from build_biomod_functions
  # Save the model and projection to personal computer
  #
  # Save out the model object
  EC_save (model_sdm, name = "model.object.RData")

  # Remove the current climate rasters to release disk space
  EC_raster_remove (predictor_info$current_climate_orig)

  # Save the projection
  EC_save_projection (model.proj, model_options_geographical$projection_name,
                      model_options_geographical$species_name,
                      model_options_geographical$species_algo_str)
}

