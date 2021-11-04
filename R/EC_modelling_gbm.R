#' Function to run an Generalized Boosting Model (GBM) within EcoCommons. It
#' is a machine learning model that uses a combo of decision trees and boosting
#' to iteratively fit random subsets of data so that each new tree takes the
#' error of the last trees into account
#'
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param mode_compute Model data parameters; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 getFormalModel
#' @importFrom biomod2 variables_importance
#' @importFrom ggplot2 ggplot
#' 
#' @export EC_modelling_gbm


EC_modelling_gbm <- function(EC.params,       # EC.params
                             response_info,   # from EC_build_response
                             model_compute) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'GBM'
  # specific parameters to run gbm algorithm
  model_options_gbm <- EC_options_gbm (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (GBM = model_options_gbm)

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
  
  EC_plot_VIP_gbm (method = model_algorithm, data1 = data1, pdf = TRUE,
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

} # end EC_modelling_gbm()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_gbm()

EC_options_gbm <- function(a){
  # Set specific parameters to run gbm algorithm
  list( distribution      = a$distribution, # "bernoulli", "gaussian", "laplace", "tdist", "huberized", "multinomial", "adaboost", "poisson", "coxph", "quantile", or "pairwise"
        n.trees           = a$n_trees, # The total number of trees to fit
        interaction.depth = a$interaction_depth, # Maximum depth of variable interactions
        n.minobsinnode    = a$n_minobsinnode,
        shrinkage         = a$shrinkage, # A shrinkage parameter applied to each tree in the expansion
        bag.fraction      = a$bag_fraction, # The fraction of the training set observations randomly selected to propose the next tree in the expansion
        train.fraction    = a$train_fraction, # The first train.fraction * nrows(data) observations are used to fit the gbm
        cv.folds          = a$cv_folds # Number of cross-validation folds to perform
  )
}


#_____________________________________________________________________________

EC_plot_VIP_gbm <- function (fittedmodel    = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
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
  
  # select the full model generated
  filekeep <-  paste(this.dir, "/", filenames[1], sep= "")
  
  working <- load(filekeep)
  fittedmodel <- biomod2::getFormalModel(eval(parse(text=working)))
  
  # variable importance plot using the inbuilt biomod2 function 'variables_importance'
  nd = dim(data1)[2]
  RespV1 = data1[,1]; subdata1 = data1[,2:nd]
  
  vi_biomod = biomod2::variables_importance(fittedmodel,data=subdata1)$mat
  nx = length(vi_biomod)
  dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
  dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
  pv <- ggplot2::ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
    labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
  ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
  
  EC_save_pdf(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
}