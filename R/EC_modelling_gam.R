#' Function to run a Generalized Additive Modes (GAM) within EcoCommons.
#' It is an statistic regression that builds multiple linear regression 
#' models by partiontioning the data and running linear regressions on 
#' each partition.
#' 
#' @param x a list (S3?) created from a json file, containing $params & $env (formerly EC.params & EC.env)
#' @param EC.params Parameters object; appears to be a nested named list. Consider converting to S3
#' @param EC.env Environment object; as for EC.params. Consider converting to S3
#'
#' @importFrom biomod2 BIOMOD_ModelingOptions
#' @importFrom biomod2 BIOMOD_Modeling
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_LoadModels
#' @export EC_modelling_gam


EC_modelling_gam <- function(a,# EC.params
                             response_info,  # from EC_build_response
                             predictor_info,  # from EC_build_predictor
                             dataset_info) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'GAM'
  # specific parameters to run GAM algorithm
  model_options_gam <- EC_options_gam (a)
  # general parameters to run biomod2 package modelling
  model_options_statregr <- EC_options_statregr (a, response_info,
                                                 model_algorithm)
  
  # Model accuracy statistics
  model_accuracy_gam <- c(model_options_statregr$biomod_eval_method,
                          model_options_statregr$dismo_eval_method)
  
  # Determine the number of pseudo absence points from pa_ratio
  pa_ratio <- EC.params$pa_ratio
  pa_number_point <- 0
  if (pa_ratio > 0) {
    pa_number_point <- floor(pa_ratio * nrow(occur))
  }
  
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions(GAM = model_options_gam)
  
  # Define model options and compute the model
  # uses biomod2 model options
  model_compute <- EC_compute_statregr (predictor_info, pa_number_point, a,
                                        model_options_statregr, model_options_gam)
  
  
  # Save a variable importance plot (VIP) and the model object
  model_save <- EC_save_statregr_model (model_compute, model_algorithm,
                                        model_options_statregr)
  
  # Project over climate scenario. Also convert projection output from grd to
  # gtiff and save the projection. Two options:
  
  # 1. projection without constraint if all env data layers are continuous
  model_projection_constrained <-
    project_model_current_constrained (predictor_info,model_compute,
                                       model_options_statregr, response_info,
                                       a, model_algorithm)
  
  # 2. projection with constraint
  model_projection <- project_model_current (model_compute, predictor_info,
                                             model_options_statregr, response_info,
                                             a, model_algorithm)

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
  