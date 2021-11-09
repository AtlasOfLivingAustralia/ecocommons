#' Function to run a Generalized Linear Model (GLM) within EcoCommons. It is an
#' statistic regression, fitted with maximum likelihood estimation for data 
#' with non-normal distribution
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
#' @importFrom tidyr pivot_longer
#' 
#' @export EC_modelling_glm


EC_modelling_glm <- function(EC.params,       # EC.params
                             response_info,   # from EC_build_response
                             model_compute) {  # from EC_build_dataset 
  
  # Set parameters to perform modelling
  model_algorithm <- 'GLM'
  # specific parameters to run GLM algorithm
  model_options_glm <- EC_options_glm (EC.params)
  # general parameters to run biomod2 package modelling
  model_options_algorithm <- EC_options_algorithm (EC.params, response_info,
                                               model_algorithm)
  # Define the model options
  model_options <- biomod2::BIOMOD_ModelingOptions (GLM = model_options_glm)
  
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
  
  EC_plot_VIP_glm (fittedmodel = model_sdm,
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
  
} # end EC_modelling_glm()



#________________________________________________________________
# subfunctions to run EC_modelling_glm()

EC_options_glm <- function(a){
  # Set specific parameters to run GLM algorithm
  list(
  	type      = a$type,	 # "simple", "quadratic" or "polynomial"; off if myFormula is not NULL
  	interaction_level = a$interaction_level,  # integer; off if myFormula is not NULL
  	# myFormula = NULL,  # specific formula; if not NULL, type and interaction.level are args are switched off
  	test      = a$test,  # "AIC", "BIC" or "none"
  	family    = a$family,  # e.g. "binomial", "gaussian", "gamma", "inverse.gaussian"
  	mustart   = a$mustart,  # starting values for the vector of means
  	control   = list(
  		epsilon = a$control_epsilon,  # positive convergence tolerance
  		maxit   = a$control_maxit,  # integer giving the maximal number of IWLS iterations
  		trace   = a$control_trace  # logical indicating if output should be produced for each iteration
  	)
  )
}


#_____________________________________________________________________________


EC_plot_VIP_glm <- function (fittedmodel  = NULL,  # or object obtained from biomod2 "BIOMOD_Modelling"
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
  
  fitteddata = fittedmodel$data  # these are the data used by biomod2 for model fitting
  nd         = dim(fitteddata)[2]
  sub.data   = fitteddata[,2:nd, drop=FALSE]
  # Can only scale numeric data
  cat.data = Filter(is.factor, sub.data)
  sub.data = Filter(is.numeric, sub.data)
  RespV    = fitteddata[,1, drop=FALSE]
  rescaled.data <- scale(sub.data)
  
  # attributes(rescaled.data)
  all.data <- cbind (RespV, as.data.frame(rescaled.data), cat.data)
  
  # head(all.data); dim(all.data)
  rescaled.glm <- glm (fittedmodel$formula, data= all.data, family= binomial)
  
  scaled = as.data.frame(cbind(coef(rescaled.glm), confint(rescaled.glm)))
  scaled = scaled[-1,]
  raw = as.data.frame(cbind(coef(fittedmodel), confint(fittedmodel)))
  raw = raw[-1,]
  
  #  variable importance plot in terms of relative effect size
  #  the relative effect size is defined as the coefficient estimated based on scaled data
  nx = length(coef(fittedmodel)[-1])
  df1 = as.data.frame(cbind(1:nx, round(scaled, 3)))
  
  names(df1) = c("xvar", "meanest", "lower", "upper")
  df1$xvar = factor(df1$xvar, labels = rownames(df1))
  
  p1 <- ggplot2::ggplot(df1, aes(x= xvar, y= meanest)) + 
    geom_hline(yintercept = 0) + 
    labs(x= " ")
  
  ps = p1 + geom_errorbar(aes(ymin= lower,ymax= upper),lwd= 0.8,width= 0.25) +
    labs(y="relative effect size") +  labs(title="         scaled data") +
    coord_flip()
  
  df2 = as.data.frame(cbind(1:nx, round(raw, 3)))
  names(df2) = c("xvar", "meanest", "lower", "upper")
  df2$xvar = factor(df2$xvar, labels = rownames(df2))
  
  if (output.table)
  {
    df1t = df1; df2t = df2
    names(df1t) = c("x.var", "coeff.est.scaled", "lower", "upper")
    names(df2t) = c("x.var", "coeff.est.raw", "lower", "upper")
    dfout = cbind (df2t,df1t)
    write.csv(dfout,file= paste(filekeep,"paraest_out.csv", sep= "_"),row.names= FALSE)
  }
  
  #  the heatmap in terms of correlation among numerical predictor variables
  rescdata = Filter(is.numeric, rescaled.glm$data[,-1, drop= FALSE])
  
  if("spearman" %in% cor.method) {
    xx = cor(rescdata, method= "spearman")
  } else if("pearson" %in% cor.method)  {
    xx = cor(rescdata)
  }
  
  lower_tri <- xx
  lower_tri[upper.tri(lower_tri)] <- NA
  
  xx.ml <- tidyr::pivot_longer(lower_tri, na.rm= TRUE)
  
  corx = xx.ml[,3]
  rm = which(is.na(corx) == TRUE)
  xx.ml = xx.ml[-rm,]
  
  pheat <- ggplot2::ggplot(xx.ml, aes(X1, X2)) + geom_tile(aes(fill = value), colour= "black") +
    scale_fill_gradient2(low = "green4", high = "violetred", mid= "white",
                         midpoint=0, limit=c(-1,1)) + labs(y=" ") + theme_minimal() +
    scale_x_discrete(limits= rownames(xx)) + scale_y_discrete(limits= colnames(xx)) + 
    coord_fixed() +
    theme(axis.title.x= element_blank(),legend.position = "bottom", axis.text.x= element_text(angle= -90)) +
    guides(fill= guide_legend(title= "correlation"))
  
  # Save as variable correlation plot.
  filename1 = sub("vip_plot", "variable_correlations", filename)
  EC_save_pdf(pheat, ncol= 1, nrow= 1, filename= filename1, aspdf= pdf)
  
  # variable importance plot in terms of AIC scores which represent the information loss,
  # e.g., the AIC score of a predictor variable representing the information loss
  # if this variable is not included in the selected model.
  
  nd = dim(data1)[2]
  
  RespV1 = data1[,1]; subdata1 = data1[,2:nd, drop=FALSE]
  glm.all = glm(formula = RespV1 ~ ., family = binomial, data = subdata1)
  
  Xaic = NULL
  for (i in 1:(nd-1))
  {
    subdf = subdata1[,-i, drop=FALSE]
    glm.one = glm(formula = RespV1 ~ . , family = binomial, data = subdf)
    Xaic = c(Xaic,AIC(glm.one))
  }
  
  relaAIC = round(Xaic - AIC(glm.all),2)
  nx = length(relaAIC)
  dfa = as.data.frame(cbind(1:nx,relaAIC))
  dfa$V1 = factor(dfa$V1, labels = rev(names(subdata1)))
  pa <- ggplot2::ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") +
    labs(y="AIC score for information loss") + labs(title="AIC approach")
  
  ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()
  
  # the variable importance plot using the inbuilt biomod2 function 'variables_importance'
  vi_biomod = biomod2::variables_importance(fittedmodel,data= subdata1)$mat
  nx = length(vi_biomod)
  dfvi = as.data.frame(cbind(1:nx, vi_biomod[,1]))
  dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
  pv <- ggplot2::ggplot(dfvi, aes(x= V1, y=rev(vi_biomod[,1]))) + labs(x= "predictor variables") +
    labs(y= "variable importance score") + labs(title= "biomod2 function 'variables_importance'")
  
  ppv = pv + geom_col(alpha= 0.6,col= "green4") + coord_flip()
  
  # Save as variable relative contribution plot.
  filename1 = sub("vip_plot", "variable_relative_contribution", filename)
  if (biom_vi) {
    EC_save_pdf(ps, ppv, ncol= 2, nrow= 1, filename= filename1, aspdf= pdf)
  }
  else {
    EC_save_pdf(ps, ppa, ncol= 2, nrow= 1, filename= filename1, aspdf= pdf)
  }
}
  
