#' subfunctions, for the Geographical algorithms scripts
#' 
#' Not exported

EC_options_geographical <- function(
  a, # formerly EC.params
  response_info, # from build_response()
  model_algorithm){
  # Build evaluation parameters for biomod2 and dismo modelling
  # 
  list(
    biomod_eval_method <- c("KAPPA", "TSS", "ROC" ,"FAR", "SR","ACCURACY","BIAS",
                            "POD", "CSI", "ETS"), #vector of evaluation metrics

    # available from dismo::evaluate.R. Not originally implemented in biomod2::Evaluate.models.R
    dismo_eval_method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR"),

    opt.scale <- ifelse(is.null(a$scale), 1000, a$scale),
    opt.tails <- a$tails, # default "both"; use to ignore the left or right
    # tail of the percentile distribution
    opt.ext <- NULL, #an optional extent object to limit the prediction to a sub-region of 'x'
    projection_name <- "current",
    species_algo_str <- if (!is.null(a$subset)) {
                              paste0("_", model_algorithm, "_",
                                     response_info$occur_species)
      }
  )
}


EC_compute_geographical <- function(
  predictor_info,  # from EC_build_predictor
  a,  # formerly EC.params
  model_options_geographical){  # from geographical_options
  # Set parameters on biomod2 format to run modelling
  
  EC_format_biomod2 (true.absen             = predictor_info$absen,
                     pseudo.absen.strategy  = a$pa_strategy,
                     pseudo.absen.disk.min  = a$pa_disk_min,
                     pseudo.absen.disk.max  = a$pa_disk_max,
                     pseudo.absen.sre.quant = a$pa_sre_quant,
                     climate.data           = predictor_info$current_climate,
                     occur                  = predictor_info$occur,
                     species.name           = model_options_geographical$species_name,
                     species_algo_str       = model_options_geographical$species_algo_str)
}





EC_save_geographical_model<-function(
  model_compute_sdm,  # from model_compute
  model_algorithm,
  model_options_geographical){  # from build_biomod_functions
  # Save the model and projection to personal computer
  #
  # Save out the model object
  EC_save(model_sdm, name = "model.object.RData")

  # Remove the current climate rasters to release disk space
  EC_raster_remove(predictor_info$current_climate_orig)

  # Save the projection
  EC_save_projection (model.proj, model_options_geographical$projection_name,
                      model_options_geographical$species_name,
                      model_options_geographical$species_algo_str)
}

