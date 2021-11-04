#' Function to run an Classification Tree Algorithm (CTA) within EcoCommons. It
#' is a machine learning model that repeatedly splits the datasets into groups
#' based on a threshold value of one of the env variables
#'
#' @param EC.params List created from a json file, containing source_file$params
#' @param response_info Response object; a nested named list created on EC_build_response
#' @param mode_compute Model data parameters; a nested named list created on EC_build_dataset
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 getFormalModel
#' @importFrom biomod2 variables_importance
#' @importFrom ggdendro dendro_data
#' @importFrom ggdendro theme_dendro
#' @importFrom ggplot2 ggplot
#' 
#' @export EC_modelling_cta


EC_modelling_cta <- function(EC.params,        # EC.params
                             response_info,    # from EC_build_response
                             model_compute) {  # from EC_build_dataset 

  # Set parameters to perform modelling
  model_algorithm <- 'CTA'
  # specific parameters to run CTA algorithm
  model_options_cta <- EC_options_cta (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)

  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (cta = model_options_cta)

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
  
  EC_plot_VIP_cta (fittedmodel = model_sdm,
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

} # end EC_modelling_cta()


#_____________________________________________________________________________
# subfunctions to run EC_modelling_cta()

EC_options_cta <- function(a){
  # Set specific parameters to run cta algorithm
  list( method = a$method, #"anova", "poisson", "class" or "exp"
        # parms = "default", #optional parameters for the splitting function
        # cost = NULL, #a vector of non-negative costs, one for each variable in the model. Defaults to one for all variables
        control = list(
          minsplit       = a$control_minsplit, #the minimum number of observations that must exist in a node in order for a split to be attempted
          minbucket      = a$control_minbucket, #the minimum number of observations in any terminal <leaf> node
          cp             = a$control_cp, #complexity parameter
          maxcompete     = a$control_maxcompete, #number of competitor splits retained in the output
          maxsurrogate   = a$control_maxsurrogate, #number of surrogate splits retained in the output
          usesurrogate   = a$control_usesurrogate, #how to use surrogates in the splitting proces
          xval           = a$control_xval, #number of cross-validations
          surrogatestyle = a$control_surrogatestyle, #controls the selection of a best surrogate
          maxdepth       = a$control_maxdepth  #Set the maximum depth of any node of the final tree, with the root node counted as depth 0. Values greater than 30 rpart will give nonsense results on 32-bit machines
          )
        )
}


#_____________________________________________________________________________

EC_plot_VIP_cta <- function (fittedmodel  = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
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
  
  # variable importance is part of the model fitting outcomes with 'rpart' algorithm and
  # this information can be used for generating the variable importance plot
  
  varimp0 = fittedmodel$variable.importance
  nx = length(varimp0)
  df0 = as.data.frame(cbind(1:nx,varimp0))
  df0$V1 = factor(df0$V1, labels = rev(names(varimp0)))
  
  p <- ggplot2::ggplot(df0, aes(x=V1, y=rev(varimp0))) + 
    labs(x=" ") +
    labs(y= "variable importance score") + 
    labs(title= "part of the 'rpart' model output")
  
  pp0 = p + geom_col(alpha=0.6,col="blue") + 
    coord_flip()
  
  ddata <- ggdendro::dendro_data(fittedmodel)
  ppt = ggplot2::ggplot() +
    geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = ddata$labels, aes(x = x, y = y, label = label), size = 3, vjust = 0) +
    geom_text(data = ddata$leaf_labels, aes(x = x, y = y, label = label), size = 3, vjust = 1) +
    ggdendro::theme_dendro()
  
  # variable importance plot using the inbuilt biomod2 function 'variables_importance'
  
  nd        = dim(data1)[2]
  subdata1  = data1[,2:nd, drop=FALSE]
  vi_biomod = biomod2::variables_importance(fittedmodel,data=subdata1)$mat
  nx        = length(vi_biomod)
  dfvi      = as.data.frame(cbind(1:nx,vi_biomod[,1]))
  dfvi$V1   = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
  pv <- ggplot2::ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) +
    labs(x="predictor variables") +
    labs(y="variable importance score") + 
    labs(title="biomod2 function 'variables_importance'")
  
  ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
  
  if (biom_vi) {
    EC_save_pdf(ppt, ppv, ncol=2, nrow=1, filename=filename, aspdf=pdf)
    
  } else {
    EC_save_pdf(ppt, pp0, ncol=2, nrow=1, filename=filename, aspdf=pdf)
  }
}