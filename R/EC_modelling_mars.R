#' Function to run a Multivariate Adaptive Regression Splines (MARS) within 
#' EcoCommons.It is an statistic regression that builds Builds multiple linear 
#' regression models by partintioning the data and running linear regressions 
#' on each partition
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
#' @export EC_modelling_mars


EC_modelling_mars <- function( EC.params,       # EC.params
                               response_info,   # from EC_build_response
                               model_compute) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'MARS'
  # specific parameters to run MARS algorithm
  model_options_mars <- EC_options_mars (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (MARS = model_options_mars)
 
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
  x.data <- attr(model_compute, "data.env.var")
  y.data <- attr(model_compute, "data.species")
  data1 <- data.frame(y.data, x.data)
  
  EC_plot_VIP_mars (fittedmodel = model_sdm,
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

} # end EC_modelling_mars()



#________________________________________________________________
# subfunctions to run EC_modelling_mars()

EC_options_mars <- function(a){
  # Set specific parameters to run GAM algorithm
  list(
    type     = a$type, #"simple", "quadratic" or "polynomial"; switched off if myFormula is not NULL
    interaction.level = a$interaction_level,   # integer; interaction between variables. Switched off if myFormula is not NULL
    # myFormula = NULL, #specific formula; if not NULL, type and interaction.level args are switched off (N/A for EcoCommons)
    nk       =  a$nk, # maximum number of model terms
    penalty  = a$penalty, #generalized cross validation (gcv) penalty per knot
    thresh   = a$thresh, #forward stepwise stopping threshold
    nprune   = a$nprune, #maximum number of terms in the pruned model
    pmethod  = a$pmethod #pruning method
  )
}



#_____________________________________________________________________________


EC_plot_VIP_mars <- function (fittedmodel  = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
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
  fittedmodel <- biomod2::BIOMOD_LoadModels(working)
  
  # variable importance plot using the inbuilt function 'varImp' from package 'caret'
  nd      = dim(data1)[2]
  RespV1  = data1[,1]; subdata1 = data1[,2:nd]
  var_imp = varImp(fittedmodel)
  nx      = length(var_imp[,1])
  dfvi    = as.data.frame(cbind(1:nx,var_imp[,1]))
  dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(var_imp)))
  pv <- ggplot2::ggplot(dfvi, aes(x= V1, y= rev(var_imp[,1]))) + 
    labs(x="predictor variables") +
    labs(y="relative reduction in GCV") + 
    labs(title="function 'varImp' in package 'caret'")
  
  ppv = pv + geom_col(alpha = 0.6, col ="red") + coord_flip()
  
  EC_save_pdf(ppv, ncol = 1, nrow = 1, filename = filename, aspdf = pdf)
}