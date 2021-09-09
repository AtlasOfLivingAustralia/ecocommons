#' Function to run an Boosted Regression Tree (BRT) within EcoCommons. It
#' is a machine learning model that uses a combo of decision trees and boosting
#' to iteratively fit random subsets of data so that each new tree takes the
#' error of the last trees into acocunt
#'
#' @param a a list created from a json file, containing $params & $env (formerly EC.params & EC.env)
# @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
# @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @export EC_modelling_brt


EC_modelling_brt <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset

  # Set parameters to perform modelling
  model_algorithm <- 'BRT'
  # specific parameters to run brt algorithm
  model_options_brt <- EC_options_brt (a)

  # Model accuracy statistics
  biomod_eval_method <- c("KAPPA", "TSS", "ROC" ,"FAR", "SR","ACCURACY","BIAS",
                          "POD", "CSI", "ETS"), #vector of evaluation metrics

  # available from dismo::evaluate.R. Not originally implemented in biomod2::Evaluate.models.R
  dismo_eval_method <- c("ODP", "TNR", "FPR", "FNR", "NPP", "MCR", "OR"),

  model_accuracy_brt <- c(biomod_eval_method, dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <-  EC_format_biomod2 (true.absen             = predictor_info$absen,
                                       pseudo.absen.strategy  = a$pa_strategy,
                                       pseudo.absen.disk.min  = a$pa_disk_min,
                                       pseudo.absen.disk.max  = a$pa_disk_max,
                                       pseudo.absen.sre.quant = a$pa_sre_quant,
                                       climate.data           = predictor_info$current_climate,
                                       occur                  = predictor_info$occur,
                                       species.name           = model_options_geographical$species_name,
                                       species_algo_str       = model_options_geographical$species_algo_str)

  # Extract occurrence and absence data
  coord = cbind(model_compute@coord, model_compute@data.env.var)
  occur = coord[c(which(model_compute@data.species == 1)), names(coord)]
  absen = coord[c(which(model_compute@data.species == 0 | is.na(model_compute@data.species))), names(coord)]


  brt.data = coord
  # setup the data as needed
  brt.data$pa = c(rep(1,nrow(occur)),rep(0,nrow(absen)))
  # run the algorithm
  model.sdm <- gbm.step(data = brt.data,
                        gbm.x = which(names(brt.data) %in% names(predictor_info$current_climate)),
                        gbm.y=which(names(brt.data)=='pa'),
                        fold.vector = brt.fold.vector,
                        tree.complexity = brt.tree.complexity,
                        learning.rate = brt.learning.rate,
                        bag.fraction = brt.bag.fraction,
                        site.weights = brt.site.weights,
                        var.monotone = brt.var.monotone,
                        n.folds = brt.n.folds,
                        prev.stratify = brt.prev.stratify,
                        family = brt.family,
                        n.trees = brt.n.trees,
                        step.size = brt.step.size,
                        max.trees = brt.max.trees,
                        tolerance.method = brt.tolerance.method,
                        tolerance = brt.tolerance,
                        keep.data = brt.keep.data,
                        plot.main = brt.plot.main,
                        plot.folds = brt.plot.folds,
                        verbose = brt.verbose,
                        silent = brt.silent,
                        keep.fold.models = brt.keep.fold.models,
                        keep.fold.vector = brt.keep.fold.vector,
                        keep.fold.fit = brt.keep.fold.fit)



  # Save the model object
  EC_save (model_sdm, name = "model.object.RData")

  # Prediction over current climate cenario

  # 1. projection without constraint if all env data layers are continuous
  if (predictor_info$genUnconstraintMap &&
      all(predictor_info$type == 'continuous') &&
      (!is.null(constraint_info$constraints) || predictor_info$generateCHull)) {
    model.proj <-
      predict(dataset_info$current_climate_orig,
              model_sdm, n.trees=model_sdm$gbm.call$best.trees, type="response")

    # remove the current_climate to release disk space
    EC_raster_remove (predictor_info$current_climate)

    # save the projection
    EC_save_projection (model.proj, projection_name,
                        response_info$occur_species, species_algo_str,
                        filename_ext="unconstrained")
  }

  # 2. projection with constraint
  model.proj = predict(predictor_info$current_climate, model_sdm,
                       n.trees=model.sdm$gbm.call$best.trees, type="response")


  #remove the current_climate to release disk space
  EC_raster_remove (predictor_info$current_climate)

  # save the projection
  EC_save_projection (model.proj, projection_name,
                      response_info$occur_species, species_algo_str)

  # evaluate model
  EC_save_dismo_eval ('BRT', model_sdm, occur, absen,
                      response_info$occur_species)

  # RETURN?? not set yet

} # end EC_modelling_brt()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_brt()

EC_options_brt <- function(a){
  # Set specific parameters to run brt algorithm
  list( brt.fold.vector = NULL #a fold vector to be read in for cross validation with offsets
        brt.tree.complexity = a$tree_complexity #sets the complexity of individual trees
        brt.learning.rate = a$learning_rate #sets the weight applied to individual trees
        brt.bag.fraction = a$bag_fraction #sets the proportion of observations used in selecting variables
        brt.site.weights = NULL # rep(1, nrow(data)) #allows varying weighting for sites
        brt.var.monotone = NULL # rep(0, length(gbm.x)) #restricts responses to individual predictors to monotone
        brt.n.folds = a$n_folds #number of folds
        brt.prev.stratify = a$prev_stratify #prevalence stratify the folds - only for presence/absence data
        brt.family = a$family #family - bernoulli (=binomial), poisson, laplace or gaussian
        brt.n.trees = a$n_trees #number of initial trees to fit
        brt.step.size = brt.n.trees #numbers of trees to add at each cycle
        brt.max.trees = a$max_trees #max number of trees to fit before stopping
        brt.tolerance.method = a$tolerance_method #method to use in deciding to stop - "fixed" or "auto"
        brt.tolerance = a$tolerance_value #tolerance value to use - if method == fixed is absolute, if auto is multiplier * total mean deviance
        brt.keep.data = FALSE #Logical. keep raw data in final model
        brt.plot.main = FALSE #Logical. plot hold-out deviance curve
        brt.plot.folds = FALSE #Logical. plot the individual folds as well
        brt.verbose = FALSE #Logical. control amount of screen reporting
        brt.silent = FALSE #Logical. to allow running with no output for simplifying model)
        brt.keep.fold.models = FALSE #Logical. keep the fold models from cross valiation
        brt.keep.fold.vector = TRUE #Logical. allows the vector defining fold membership to be kept
        brt.keep.fold.fit = FALSE #Logical. allows the predicted values for observations from cross-validation to be kept
        projection_name = "current",
        species_algo_str = ifelse(is.null(a$subset),
                                  sprintf(model_algorithm,"%s_",
                                          response_info$occur_species),
                                  sprintf(model_algorithm,"%_%s",
                                          response_info$occur_species, a$subset))
        )
  )
}
