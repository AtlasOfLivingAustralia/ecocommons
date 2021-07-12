# Define helper functions to print parameters
#
# The EcoCommons plattform computes many Species Distribution Modelling
# algorithms to be selected by the user. Each algorithm requires a specific
# set of parameters. This funtion simply aims to compute simple names for
# parameters and assign the necessary parameters to run specific algorithm.

# Check if libraries are installed, install if necessary and then load them
necessary <- c("ggplot2","tools", "rjson", "dismo", "gbm", "rgdal", "rgeos", "pROC", "png", "gstat", "gdalUtils", "ggdendro", "raster","biomod2","rasterVis") #list the libraries needed

###REMOVED NOW JUST FOR TEST:::: "SDMTools",‘spatial.tools’ and ‘rmaxent’

installed <- necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >= 1) {
  install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
  library(lib,character.only = T) #load the libraries
}

# Load the parameters
# The 'param.json' file used here is generated on EcoCommons python script
params <- rjson::fromJSON(file="~/Documents/BCCVL_scripts/test_script.json")
EC.params <- params$params
EC.env <- params$env
rm(params)





EC_RenameParameters <- function (param, value) {
  # Computes simple names for different parameters
  #
  # Args:
  #  param: list of possible parameters that can be implemented
  #  value: ?
  #
  # Return:
  #  pname: new names for the different parameters, including value
  pname <- gsub("_", " ", param)
  if (param == "prevalence") {
    pname = "weighted response weights"
  }
  else if (param == "var_import") {
    pname = "resampling"
  }
  else if (param == "nbcv") {
    pname = "NbCV"
  }
  else if (param == "n_trees") {
    pname = "trees added each cycle"
  }
  else if (param == "control_xval") {
    pname = "cross-validations"
  }
  else if (param == "control_minbucket") {
    pname = "minimum bucket"
  }
  else if (param == "control_minsplit") {
    pname = "minimum split"
  }
  else if (param == "control_cp") {
    pname = "complexity parameter"
  }
  else if (param == "control_maxdepth") {
    pname = "maximum depth"
  }
  else if (param == "irls_reg") {
    pname = "irls.reg"
  }
  else if (param == "maxit") {
    pname = "maximum iterations"
  }
  else if (param == "mgcv_tol") {
    pname = "convergence tolerance"
  }
  else if (param == "mgcv_half") {
    pname = "number of halvings"
  }
  else if (param == "n_minobsinnode") {
    pname = "Min observations in terminal node"
  }
  else if (param == "control_epsilon") {
    pname = "control: epsilon"
  }
  else if (param == "control_maxit") {
    pname = "control: maxit"
  }
  else if (param == "control_trace") {
    pname = "control: trace"
  }
  else if (param == "model") {
    pname = "Model returned"
  }
  else if (param == "x") {
    pname = "x returned"
  }
  else if (param == "y") {
    pname = "y returned"
  }
  else if (param == "qr") {
    pname = "QR returned"
  }
  else if (param == "singular_ok") {
    pname = "Singular fit ok"
  }
  else if (param == "thresh") {
    pname = "threshold"
  }
  else if (param == "maximumiterations") {
    pname = "Maximum iterations"
  }
  else if (param == "ntree") {
    pname = "number of trees"
  }
  else if (param == "mtry") {
    pname = "number of variables at each split (mtry)"
  }
  else if (param == "nodesize") {
    pname = "node size"
  }
  else if (param == "maxnodes") {
    pname = "maximum nodes"
  }
  else if (param == "pa_ratio") {
    pname = "absence-presence ratio"
  }
  return(paste(pname, " = ", value, "\n", sep="", collapse=""))
}

EC_ParameterPrint <- function(params) {
  # Assign parameters necessary to run specific algorithms
  #
  # Args:
  #  parameters: set of parameters that can be used
  #  func: functions or algorithms that the parameters are assigned to
  #
  # Return:
  #  list of parameters restricted to specific algorithms
  func <- params[["function"]]
  if (is.null(func))
    return("")
  cat("Algorithm:", func, "\n")

  pnames = c("random_seed")
  if (func == "ann") {
    pnames = c("prevalence", "var_import", "maxit", "nbcv", "rang", "random_seed")
  }
  else if (func == "brt") {
    pnames = c("tree_complexity", "learning_rate", "bag_fraction", "n_folds", "prev_stratify", "family", "n_trees", "max_trees", "tolerance_method", "tolerance_value", "random_seed")
  }
  else if (func == "cta") {
    pnames = c("prevalence", "var_import", "method", "control_xval", "control_minbucket", "control_minsplit", "control_cp", "control_maxdepth", "random_seed")
  }
  else if (func == "fda") {
    pnames = c("prevalence", "var_import", "method", "random_seed")
  }
  else if (func == "gam") {
    pnames = c("prevalence", "var_import", "interaction_level", "family", "irls_reg", "epsilon", "maxit", "mgcv_tol", "mgcv_half", "random_seed")
  }
  else if (func == "gamlss") {
    pnames = c("sigma_formula", "nu_formula", "tau_formula", "family", "weights", "contrasts", "method", "start_from", "mu_start", "sigma_start", "nu_start", "tau_start", "mu_fix", "sigma_fix", "nu_fix", "tau_fix", "control", "i_control", "other_args", "random_seed")
  }
  else if (func == "gbm") {
    pnames = c("prevalence", "var_import", "distribution", "n_trees", "interaction_depth", "n_minobsinnode", "shrinkage", "bag_fraction", "train_fraction", "cv_folds", "random_seed")
  }
  else if (func == "glm") {
    pnames = c("prevalence", "var_import", "type", "interaction_level", "test", "family", "mustart", "control_epsilon", "control_maxit", "control_trace", "random_seed")
  }
  else if (func == "lm") {
    pnames = c("subset", "weights", "na_action", "method", "model", "x", "y", "qr", "singular_ok", "contrasts", "offset", "random_seed")
  }
  else if (func == "manova") {
    pnames = c("projections_returned", "qr", "contrasts", "subset", "weights", "na_action", "random_seed")
  }
  else if (func == "mars") {
    pnames = c("prevalence", "var_import", "degree", "nk", "penalty", "thresh", "prune", "random_seed")
  }
  else if (func == "maxent") {
    pnames = c("prevalence", "var_import", "maximumiterations", "linear", "quadratic", "product", "threshold", "hinge", "lq2lqptthreshold", "lq2lqthreshold", "hingethreshold", "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "defaultprevalence", "random_seed")
  }
  else if (func == "rf") {
    pnames = c("prevalence", "var_import", "do.classif", "ntree", "mtry", "nodesize", "maxnodes", "random_seed")
  }
  else if (func == "sre") {
    pnames = c("prevalence", "var_import", "quant", "random_seed")
  }

  pnames = c(pnames, "pa_ratio", "pa_strategy", "pa_sre_quant", "pa_disk_min", "pa_disk_max")
  for (p in pnames) {
    cat(EC_RenameParameters(p, params[[p]]))
  }
  return("")
}


# Print out parameters used
EC_ParameterPrint(EC.params)

# print a warning when there is a problem logging the data
EC_LogWarning <-function(str, prefix="EcoCommons Warning: ") {
  print(paste(prefix, str, sep=""))
}


