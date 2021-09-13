#' Function to run BioClim within EcoCommons. It is a
#' profile model to define a multidimensional space bounded by the min and max
#' values of environment variables for all occurrences as the species potential
#' range.
#'
#' @param a formerly EC.params
#' @param response_info from EC_build_response
#' @param predictor_info from EC_build_predictor
#' @param dataset_info from EC_build_dataset
#'
#' @importFrom dismo bioclim
#' @importFrom raster writeRaster
#'
#' @export EC_modelling_bioclim

EC_modelling_bioclim <- function( a,  # EC.params
                                  response_info,  # from EC_build_response
                                  predictor_info,  # from EC_build_predictor
                                  dataset_info) {  # from EC_build_dataset

  # Set parameters to perform modelling
  model_algorithm <- 'bioclim'

  opt.tails = a$tails  # to ignore left/right tail of the percentile distribution
  opt.ext = NULL  # optional extent object to limit the prediction to a sub-region

  projection_name = "current"
  species_algo_str = ifelse(is.null(a$subset),
                            sprintf(model_algorithm,"%s_",
                                    response_info$occur_species),
                            sprintf(model_algorithm,"%_%s",
                                    response_info$occur_species, a$subset))

  biomod_eval_method <- c("KAPPA", "TSS", "ROC" ,"FAR", "SR", "ACCURACY",
                          "BIAS", "POD", "CSI", "ETS") #vector of evaluation metrics

  # available from dismo::evaluate.R. Not originally implemented in biomod2::Evaluate.models.R
  dismo_eval_method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR")

  # Model accuracy statistics
  model_accuracy_bioclim <- c(biomod_eval_method, dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio = a$pa_ratio
  pa_number_point = 0
  if (pa_ratio > 0) {
    pa_number_point = floor(pa_ratio * nrow(occur))
  }

  # 1. Format the data as required by the biomod package
  model_data <- EC_format_biomod2 ( true.absen             = predictor_info$absen,
                                    pseudo.absen.points    = pa_number_point,
                                    pseudo.absen.strategy  = a$pa_strategy,
                                    pseudo.absen.disk.min  = a$pa_disk_min,
                                    pseudo.absen.disk.max  = a$pa_disk_max,
                                    pseudo.absen.sre.quant = a$pa_sre_quant,
                                    climate.data           = predictor_info$current_climate,
                                    occur                  = predictor_info$occur,
                                    species.name           = a$species_name,
                                    species_algo_str       = species_algo_str)

  # Extract occurrence and absence data
  coord = cbind(model_data@coord, model_data@data.env.var)
  occur = coord[c(which(model_data@data.species == 1)), names(coord)]
  absen = coord[c(which(model_data@data.species == 0 | is.na(model_data@data.species))),
                names(coord)]

  # BioClim specific requirements

  if (!all(enviro.data.type=="continuous")) {
    stop("bioclim not run because categorical data cannot be used")
  } else {
    # run bioclim with matrix of enviro data
    model_sdm = bioclim(x = occur[,names(predictor_info$current_climate), drop=FALSE])
    # save out the model object
    EC_save (model_sdm, name = "model.object.RData")

    # Do projection over current climate scenario without constraint
    if (constraint_info$genUnconstraintMap &&
        (!is.null(constraint_info$constraints) || constraint_info$generateCHull)) {

      model_proj = predict(model_sdm, current.climate.scenario.orig,
                           tails=opt.tails)


      # Remove the current climate rasters to release disk space
      EC_raster_remove(predictor_info$current_climate_orig)

      # Save the projection
      EC_save_projection (model_proj, projection_name,
                          response_info$occur_species, species_algo_str,
                          filename_ext = "unconstrained")
    }


      # Predict for given climate scenario
      model_proj <- predict (model_sdm, predictor_info$current_climate[[1]],
                            tails = opt.tails)

      # remove the current climate rasters to release disk space
      EC_raster_remove (predictor_info$current_climate_orig)

      # save output
      EC_save_model_proj (model_proj, projection_name,
                          response_info$occur_species, species_algo_str)

    # evaluate model
    if (!is.null(absen)) {
      EC_save_dismo_eval ('bioclim', model_sdm, occur, absen,
                                     response_info$occur_species)
    }
  }
}


#_____________________________________________________________________________
## Subfunction, as appear within EC_modelling_bioclim()


EC_save_model_proj <- function(model_obj,
                               projection_name,
                               species,
                               species_algo_str,
                               outputdir = EC.env$outputdir, 
                               filename_ext = NULL) {
  
  filename = EC_file_path("proj", projection_name, species_algo_str, outputdir, 
                          filename_ext, "tif")
  
  raster::writeRaster (model_obj, filename, format="GTiff",
                       options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  
  pngfilename = EC_file_path("proj", projection_name, species_algo_str, 
                             outputdir, filename_ext, "png")
  
  png(pngfilename)
  
  title = paste(species, projection_name, "projections", sep=" ")
  
  plot(model_obj, xlab="latitude", ylab="longtitude", main=title)
  
  dev.off()
  
  return (filename)
}

