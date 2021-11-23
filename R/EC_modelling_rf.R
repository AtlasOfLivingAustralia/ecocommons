#' Function to run a Random Forest (RF) within EcoCommons. It is a machine
#' learning model that grows numerous decision trees on random subsets and
#' averages the predictions of the trees to estimate the importance of each
#' environment variable
#'
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param mode_compute Model data parameters; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 getFormalModel
#' @importFrom ggplot2 ggplot
#' 
#' @export EC_modelling_rf


EC_modelling_rf <- function(EC.params,       # EC.params
                            response_info,   # from EC_build_response
                            model_compute) {  # from EC_build_dataset 

  # Set parameters to perform modelling
  model_algorithm <- 'RF'
  # specific parameters to run rf algorithm
  model_options_rf <- EC_options_rf (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (RF = model_options_rf)

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
                              rescal.all.models  = model_options_algorithm$rescale_all_models,
                              do.full.models     = model_options_algorithm$do_full_models,
                              modeling.id        = model_options_algorithm$model_id)
  
  # save out the model object
  EC_save (model_sdm, name = "model.object.RData")
  
  
  # Generate and save a variable importance plot (VIP) and the model object
  x.data <- attr(model_compute, "data.env.var")
  y.data <- attr(model_compute, "data.species")
  data1 <- data.frame(y.data, x.data)
  
  EC_plot_VIP_rf (fittedmodel = model_sdm,
                  data1 = data1,
                  pdf = TRUE,
                  filename = paste('vip_plot', model_options_algorithm$species_algo_str,
                                   sep = "_"),
                  this.dir = EC.env$outputdir)
  
  
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

} # end EC_modelling_rf()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_rf()

EC_options_rf <- function(a){
  # Set specific parameters to run rf algorithm
  list( do.classif =  a$do.classif,
        ntree      =  a$ntree,
        mtry       = if (a$mtry == "default") a$mtry else as.integer(a$mtry),
        nodesize   = a$nodesize,
        maxnodes   = a$maxnodes
  )
}


#_____________________________________________________________________________

EC_plot_VIP_rf <- function (fittedmodel  = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
                            cor.method   = c("pearson", "spearman"),  # 'person' for linear data; 'spearman' for non-linear; rank-based
                            pdf          = TRUE,
                            biom_vi      = FALSE,  # function/algorithm other than biomod2 "variables_importance"
                            output.table = FALSE,  # csv file with GLM parameters and the 95% confidence interval
                            data1,  # data frame with response and predictors variables; should contain all predict variables
                            this.dir,  # route to access biomod2 model (not the model name)
                            filename) {  # to be saved without the file extension
  
  data1$y.data[is.na(data1$y.data)] <- 0
  
  # select the full model generated
  filekeep <-  paste(this.dir, "model.object.RData", sep= "/")
  
  working <- load(filekeep)
  fittedmodel <- biomod2::BIOMOD_LoadModels(model_sdm)
  
  # Random forests (rf) provide an improvement over bagged trees by way of a small tweak
  #  that decorrelates the trees.
  # Note that the variable importance plot using the inbuilt biomod2 function 'variables_importance'
  #  does not seem working with Random forests algorithm.  On the other hand, however, the fitted
  #  rf model object contains the variable importance information which is measured by the mean derease
  #  in Gini index (expressed relative to the maximum).
  # While RSS is used for measuring the regression tree model performance, the Gini index is used for
  #  measuring the classification tree model performance and Gini index is a measure of total variance
  #  across the K classes.
  
  nd      = dim(data1)[2]
  RespV1  = data1[,1]; subdata1 = data1[,2:nd]
  out.rf  = fittedmodel$importance
  rfImp   = out.rf[,1]
  nx      = length(rfImp)
  dfrf    = as.data.frame(cbind(1:nx,rfImp))
  dfrf$V1 = factor(dfrf$V1, labels = rev(names(rfImp)))
  
  prf <- ggplot2::ggplot(dfrf, aes(x=V1, y=rev(rfImp))) + 
    labs(x="predictor variables") +
    labs(y="mean decrease in Gini index") + 
    labs(title="part of rf model fitting outputs")
  
  pprf = prf + geom_col(alpha = 0.6, col = "blue") + coord_flip()
  
  EC_save_pdf(pprf, ncol = 1, nrow = 1, filename = filename, aspdf = pdf)
}