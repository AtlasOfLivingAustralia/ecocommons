#' Function to run an Classification Tree Algorithm (CTA) within EcoCommons. It
#' is a machine learning model that repeatedly splits the datasets into groups
#' based on a threshold value of one of the env variables
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
#' @export EC_modelling_cta


EC_modelling_cta <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset

  # Set parameters to perform modelling
  model_algorithm <- 'CTA'
  # specific parameters to run CTA algorithm
  model_options_cta <- EC_options_cta (a)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (a, response_info,
                                                 model_algorithm)

  # Model accuracy statistics
  model_accuracy_cta <- c(model_options_algorithm$biomod_eval_method,
                          model_options_algorithm$dismo_eval_method)

  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- a$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (cta = model_options_cta)

  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_algorithm (predictor_info, pa_number_point, a,
                                        model_options_algorithm, model_options_cta)


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

} # end EC_modelling_cta()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_cta()

EC_options_cta <- function(a){
  # Set specific parameters to run cta algorithm
  list( method = a$method, #"anova", "poisson", "class" or "exp"
        # parms = "default", #optional parameters for the splitting function
        # cost = NULL, #a vector of non-negative costs, one for each variable in the model. Defaults to one for all variables
        control = list(
          minsplit = a$control_minsplit, #the minimum number of observations that must exist in a node in order for a split to be attempted
          minbucket = a$control_minbucket, #the minimum number of observations in any terminal <leaf> node
          cp = a$control_cp, #complexity parameter
          maxcompete = a$control_maxcompete, #number of competitor splits retained in the output
          maxsurrogate = a$control_maxsurrogate, #number of surrogate splits retained in the output
          usesurrogate = a$control_usesurrogate, #how to use surrogates in the splitting proces
          xval = a$control_xval, #number of cross-validations
          surrogatestyle = a$control_surrogatestyle, #controls the selection of a best surrogate
          maxdepth = a$control_maxdepth  #Set the maximum depth of any node of the final tree, with the root node counted as depth 0. Values greater than 30 rpart will give nonsense results on 32-bit machines
          )
        )
}
