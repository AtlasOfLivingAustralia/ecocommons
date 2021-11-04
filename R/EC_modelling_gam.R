#' Function to run a Generalized Additive Modes (GAM) within EcoCommons.
#' It is an statistic regression that builds multiple linear regression 
#' models by partiontioning the data and running linear regressions on 
#' each partition.
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
#' @importFrom mgcv gam
#' 
#' @export EC_modelling_gam


EC_modelling_gam <- function(EC.params,       # EC.params
                             response_info,   # from EC_build_response
                             model_compute) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'GAM'
  # specific parameters to run GAM algorithm
  model_options_gam <- EC_options_gam (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                                 model_algorithm)
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions(GAM = model_options_gam)
  
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
  
  EC_plot_VIP_gam (method = model_algorithm, data1 = data1, pdf = TRUE,
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

} # end EC_modelling_gam()



#________________________________________________________________
# subfunctions to run EC_modelling_gam()

EC_options_gam <- function (a){
  # Set specific parameters to run GAM algorithm
  list(
    algo = "GAM_mgcv",  # options are "GAM_mgcv", "GAM_gam" or "BAM_mgcv"
    type = "s_smoother",  # used to generate the formula
    interaction.level = a$interaction_level,  # integer; interaction between variables
    myFormula = NULL, # specific formula; if not NULL, type and interaction.level are args are switched off
    k = a$k, # a smooth term in a formula argument to gam (see gam s or mgcv s)
    family = a$family, # e.g. "binomial", "bernoulli", "gaussian", "laplace"
    method = a$method, # the smoothing parameter estimation method
    
    # Set a numerical optimization method to smoothing parameter estimation 
    # criterion (given by method). "perf" (deprecated) for performance iteration;
    # "outer" for the more stable direct approach. NOTE: "outer" requires
    # alternative optimizers, specified in the second element: "newton" (default),
    # "bfgs", "optim", "nlm" and "nlm.fd"
    optimizer = a$optimizer,
    select = a$select, # If TRUE- GAM add an extra penalty to each term
    knots = NULL,  # optional; list to be used for basis construction
    paraPen = NULL,  # optional; list with penalties to be applied to parametric model terms   
    control = list(
      # Specific GAM settings
      irls.reg = a$control_irls.reg, # the size of the ridge regression penalty to the model to impose identifiability; for most models this should be 0
      epsilon = a$control_epsilon, # this is used for judging conversion of the GLM IRLS loop
      maxit = a$control_maxit, # maximum number of IRLS iterations to perform
      mgcv.tol = a$control_mgcv.tol, # the convergence tolerance parameter to use in GCV/UBRE optimization
      mgcv.half = a$control_mgcv.half, # if a step of the GCV/UBRE optimization method leads to a worse GCV/UBRE score, then the step length is halved; this is the number of halvings to try before giving up
      trace = a$control_trace, # set this to TRUE to turn on diagnostic output
      rank.tol = .Machine$double.eps^0.5, # the tolerance used to estimate the rank of the fitting problem
      nlm = list(), # list of control parameters to pass to nlm if this is used for outer estimation of smoothing parameters (not default)
      optim = list(), # list of control parameters to pass to optim if this is used for outer estimation of smoothing parameters (not default)
      newton = list(), # list of control parameters to pass to default Newton optimizer used for outer estimation of log smoothing parameters
      outerPIsteps = 0, # the number of performance interation steps used to initialize outer iteration
      idLinksBases = TRUE, # if smooth terms have their smoothing parameters linked via the id mechanism (see s), should they also have the same bases. Set this to FALSE only if you are sure you know what you are doing
      scalePenalty = TRUE, # this option rescales the penalty matrices to accomodate this problem. Probably should be set to FALSE if you are linking smoothing parameters but have set idLinkBases to FALSE
      keepData = FALSE # should a copy of the original data argument be kept in the gam object
    )
  )
}


#_____________________________________________________________________________


EC_plot_VIP_gam <- function (fittedmodel  = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
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
  
  if (biom_vi)
  {
    # variable importance plot using the inbuilt biomod2 function 'variables_importance'
    nd = dim(data1)[2]
    RespV1 = data1[,1]; subdata1 = data1[,2:nd, drop=FALSE]
    
    vi_biomod = biomod2::variables_importance(fittedmodel,data=subdata1)$mat
    nx = length(vi_biomod)
    dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
    pv <- ggplot2::ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
    ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()
    
    EC_save_pdf(ppV, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  } else {
    # variable importance plot following the AIC approach
    nd = dim(data1)[2]
    RespV1 = data1[,1]
    subdata1 = data1[,2:nd, drop=FALSE]
    
    # gam function cannot take categorical data, so exclude categorical data.
    subdata1 = Filter(is.numeric, subdata1)
    
    xname = names(subdata1)
    sname = paste("s(", xname, ")",sep="")
    
    gamformu.all <- as.formula(paste("RespV1 ~ 1 +", paste(sname, collapse= "+")))
    gam.all = mgcv::gam(formula = gamformu.all, family = binomial, data = subdata1)
    
    Xaic = NULL
    nd = dim(subdata1)[2]
    for (i in 1:nd)
    {
      subdf = subdata1[, -i, drop=FALSE]
      xname1 = names(subdf)
      sname1 = paste("s(", xname1, ")",sep="")
      gamformu1 <- as.formula(paste("RespV1 ~ 1 +", paste(sname1, collapse= "+")))
      gam.one = mgcv::gam(formula = gamformu1, family = binomial, data = subdf)
      Xaic = c(Xaic,AIC(gam.one))
    }
    
    relaAIC = round(Xaic - AIC(gam.all),2)
    nx = length(relaAIC)
    dfa = as.data.frame(cbind(1:nx,relaAIC))
    dfa$V1 = factor(dfa$V1, labels = rev(names(subdata1)))
    pa <- ggplot2::ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") +
      labs(y="AIC score for information loss") + labs(title="AIC approach")
    
    ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()
    
    EC_save_pdf(ppa, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  }
}