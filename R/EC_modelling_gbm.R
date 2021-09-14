#' Function to run an Generalized Boosting Model (GBM) within EcoCommons. It
#' is a machine learning model that uses a combo of decision trees and boosting
#' to iteratively fit random subsets of data so that each new tree takes the
#' error of the last trees into account
#'
#' @param a List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param predictor_info Predictor object; a nested named list created on EC_build_predictor
#' @param dataset_info Dataset object; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' 
#' @export EC_modelling_gbm


EC_modelling_gbm <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset

  # Set parameters to perform modelling
  model_algorithm <- 'GBM'
  # specific parameters to run gbm algorithm
  model_options_gbm <- EC_options_gbm (a)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (a, response_info,
                                                 model_algorithm)

  # Model accuracy statistics
  model_accuracy_gbm <- c(model_options_algorithm$biomod_eval_method,
                          model_options_algorithm$dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (GBM = model_options_gbm)

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_algorithm (predictor_info, pa_number_point, a,
                                        model_options_algorithm, model_options_gbm)


  # Save a variable importance plot (VIP) and the model object
  model_save <- EC_save_algorithm_model (model_compute, model_algorithm,
                                        model_options_algorithm)

  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:

  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained (predictor_info,model_compute,
                                       model_options_algorithm, response_info,
                                       a, model_algorithm)

  # 2. projection with constraint
  model_projection <- project_model_current (model_compute, predictor_info,
                                             model_options_algorithm, response_info,
                                             a, model_algorithm)

  # RETURN?? not set yet

} # end EC_modelling_gbm()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_gbm()

EC_options_gbm <- function(a){
  # Set specific parameters to run gbm algorithm
  list( distribution = a$distribution, # "bernoulli", "gaussian", "laplace", "tdist", "huberized", "multinomial", "adaboost", "poisson", "coxph", "quantile", or "pairwise"
        n.trees = a$n_trees, # The total number of trees to fit
        interaction.depth = a$interaction_depth, # Maximum depth of variable interactions
        n.minobsinnode = a$n_minobsinnode,
        shrinkage = a$shrinkage, # A shrinkage parameter applied to each tree in the expansion
        bag.fraction = a$bag_fraction, # The fraction of the training set observations randomly selected to propose the next tree in the expansion
        train.fraction = a$train_fraction, # The first train.fraction * nrows(data) observations are used to fit the gbm
        cv.folds = a$cv_folds # Number of cross-validation folds to perform
  )
}
