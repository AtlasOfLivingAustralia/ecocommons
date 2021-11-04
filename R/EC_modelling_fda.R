#' Function to run a Flexible Discriminant Analysis (FDA) within EcoCommons.
#' It is an statistic regression, a classification model based on a mixture
#' of linear regression models, which uses optimal scoring to transform the
#' response variable so that the data are in a better form for linear separation,
#' and multiple adaptive regression splines to generate the discriminant surface
#' 
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param mode_compute Model data parameters; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom ggplot2 ggplot
#' 
#' @export modelling_fda


modelling_fda <- function(EC.params,       # EC.params
                          response_info,   # from EC_build_response
                          model_compute) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'FDA'
  # specific parameters to run FDA algorithm
  model_options_fda <- EC_options_fda (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)
 
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (FDA = model_options_fda)
  
  # Use define model options (created on model_compute) and compute the model
  # uses biomod2 model options
  model_sdm <-
    biomod2::BIOMOD_Modeling (data               = model_compute,
                              models             = model_algorithm,
                              models.options     = model_options,
                              NbRunEval          = model_options_algorithm$NbRunEval,
                              DataSplit          = model_options_algorithm$DataSplit,
                              Yweights           = model_options_algorithm$Yweights,
                              Prevalence         = model_options_algorithm$Prevalence,
                              VarImport          = model_options_algorithm$VarImport,
                              models.eval.meth   = model_options_algorithm$biomod_eval_meth,
                              SaveObj            = TRUE,
                              rescal.all.models  = model_options_algorithm$rescal_all_models,
                              do.full.models     = model_options_algorithm$do_full_models,
                              modeling.id        = model_options_algorithm$model_id)
  
  # save out the model object
  EC_save (model_sdm, name = "model.object.RData")
  
  
  # Generate and save a variable importance plot (VIP) and the model object
  x.data <- attr(model_compute$model_data, "data.env.var")
  y.data <- attr(model_compute$model_data, "data.species")
  data1 <- data.frame(y.data, x.data)
  
  EC_plot_VIP_fda (method = model_algorithm, data1 = data1, pdf = TRUE,
                   filename = paste('vip_plot', EC_options_algorithm$species_algo_str,
                                    sep = "_"),
                   this.dir = paste(EC_options_algorithm$species_name,
                                    "/models/EcoCommons", sep = ""))
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained (constraint_info, predictor_info, model_sdm,
                                       dataset_info, model_options_algorithm,
                                       model_algorithm, EC.params)
  
  # 2. projection with constraint
  model_projection <- project_model_current (model_sdm, dataset_info,
                                             model_options_algorithm, model_algorithm,
                                             EC.params)

  # RETURN?? not set yet
  
} # end modelling_fda()



#________________________________________________________________
# subfunctions to run modelling_fda()

EC_options_fda <- function(a){
  # Set specific parameters to run FDA algorithm
  list(
    model.options.fda <- list(
      method = a$method  # regression method used in optimal scaling
    )
  )
}




#_____________________________________________________________________________


EC_plot_VIP_fda <- function (fittedmodel    = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
                             method       = model_sdm,
                             cor.method   = c("pearson", "spearman"),  # 'person' for linear data; 'spearman' for non-linear; rank-based
                             pdf          = TRUE,
                             biom_vi      = FALSE,  # function/algorithm other than biomod2 "variables_importance"
                             output.table = FALSE,  # csv file with GLM parameters and the 95% confidence interval
                             data1 ) {  # data frame with response and predictors variables; should contain all predict variables
  
  data1$y.data[is.na(data1$y.data)] <- 0
  
  # extract the root of filenames used by biomod2 to save model results
  filenames <- dir(paste(EC_options_algorithm$species_name,
                         "/models/EcoCommons", sep = "")) 
}
