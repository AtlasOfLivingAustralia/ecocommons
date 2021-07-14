#############################################################
###     EcoCommons script to Species Traits Analysis      ###
#############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of EcoCommons etc.
## Date : ??
## Script and data info:
##  This script runs the maximum entropy modelling (MaxEnt) algorithm 
##  (Philips et al. 2004) for Species Distribution Modelling, using the 
##  biomod2 package on R and the input data generated on the EcoCommons 
##  platform.
##
## IMPORTANT: please run the 'EcoCommons_source' script before run each
##            algorithm

library("EcoCommons")

# Load the environmental raster layers
environ.rasterstack <- EC_EnviroStack(enviro.data.current, enviro.data.type,
                                      enviro.data.layer, "highest")

# if no species, then run across all species
if (!is.null(trait.species)) {
  trait.data <- subset(trait.data, species==trait.species)
} else {
  # use the trait filename as species name for the plots' filenames
  trait.species <- file_path_sans_ext(basename(trait.data.filename))
}


# Geographically constrained modelling and merge the environmental data into trait.data
if (!is.null(trait.data)) {
    merged.result <- bccvl.trait.constraint.merge(trait.data, trait.data.params,
                                                  environ.rasterstack,
                                                  enviro.data.constraints,
                                                  enviro.data.generateCHull,
                                                  generateGeoconstraint=FALSE)
    trait.data <- merged.result$data
    trait.data.params <- merged.result$params
}

#===================================================================
## CTA
##

## Runs a Classification Tree Analysis (for categorical trait data) or a 
## Regression Tree Analysis (for continuous trait data) to test the effect 
## of selected environmental variables on species traits

library("rpart")

# Run models
# Generate a formula for each trait
formulae_cta <- bccvl.trait.gen_formulae(trait.data.params, include_rf=TRUE)
for (formula in formulae) {
    trait_name <- gsub("[_ ]", "-", trimws(formula$trait))
    trait.cta.options <- list(formula = formula(formula$formula),  # should be: trait ~ env1 + env2 + env3 etc 
                              method = ifelse(formula$type == 'continuous',
                                              'anova', 'class'),  # "class"-categorical trait data; "anova"- continuous
                              na.action = na.rpart,  # deletes info on which trait value is missing
                              model = FALSE,
                              x = FALSE,
                              y = FALSE,
                              control = list(minsplit = EC.params$control_minsplit,  # minimum number of observations in a node for a split
                                             minbucket = EC.params$control_minbucket,  # minimum number of observations in any terminal node
                                             cp = EC.params$control_cp,  # complexity parameter
                                             maxcompete = EC.params$control_maxcompete,  # number of competitor splits retained in the output
                                             maxsurrogate = EC.params$control_maxsurrogate,  # number of surrogate splits retained in the output
                                             usesurrogate = EC.params$control_usesurrogate,  # how to use surrogates in splitting process
                                             surrogatestyle = EC.params$control_surstyle,  # controls the selection of a best surrogate
                                             xval = EC.params$control_xval,  # number of cross-validations
                                             maxdepth = EC.params$control_maxdepth  # maximum depth of any node of the final tree.root node=depth 0
                                             )
                              )

    trait.cta = rpart(formula = trait.cta.options$formula,
                      data = trait.data, # data frame containing trait and env data
                      method = trait.cta.options$method,
                      na.action = trait.cta.options$na.action,
                      model = trait.cta.options$model,
                      x = trait.cta.options$x,
                      y = trait.cta.options$y,
                      control = trait.cta.options$control)

    ### Save the results as text to file for each trait
    s <- summary(trait.cta) # saved and displayed as text
    bccvl.write.text(s, paste0(trait_name, "_cta_results.txt"))
    p <- printcp(trait.cta) # saved and displayed as text
    bccvl.write.text(p, paste0(trait_name, "_cta_results.txt"), append=TRUE)

    # save the plot as png image
    png(file.path(EC.env$outputdir, paste0(trait_name, "_cta_plot.png")))
    plot(trait.cta)
    text(trait.cta)
    dev.off()


    # Save the model
    EC_Save(trait.cta, paste0(trait_name, "_cta_model.object.RData"))
}

#===================================================================
## GAM
##

## DEFINITION

library("gam")  # or mgvc? 

na_action <- get0(EC.params$na_action)
if (is.null(na_action)) {
  na_action = get("na.fail")
}

# Generate a formula for each trait
formulae_gam <- bccvl.trait.gen_formulae(trait.data.params)
for (formula in formulae) {
  trait_name <- gsub("[_ ]", "-", trimws(formula$trait))
  if (formula$type == 'ordinal') {
    fam = ocat(R=nlevels(trait.data[[trait_name]]))
  }
  else if (formula$type == 'nominal')
  {
    # env variables for the formula
    envvar <- attr(terms(formula(formula$formula)), 'term.labels')
    fam <- multinom(K=length(envvar))
  }
  else {
    fam <- EC_FamilyFromString(EC.params$family)
  }
  
  # Run the model for each trait separately
  gam.result <- gam(formula=formula(formula$formula),
                   data=trait.data,
                   family=fam,
                   weights=NULL,
                   na.action=na_action,
                   start=NULL,
                   etastart=NULL,
                   mustart=NULL,
                   method=EC.params$method,
                   model=TRUE,
                   x=FALSE,
                   y=FALSE)
  
  # Save the model to file
  EC.Save(gam.result, paste0(trait_name, "_gam_model.object.RData"))
  
  # Save result summary to a text file
  s <- summary(gam.result)
  bccvl.write.text(s, paste0(trait_name, "_gam_result_summary.txt"))
  
  # save the plot as png image
  ofilename <- paste0(trait_name, "_gam_plotgam")
  bccvl.write.image(gam.result, ofilename, "plot.gam")
}


#===================================================================
## GLM
##
# Load the library
library("MASS")
library("nnet")


## DEFINITION

projection.name <- "current"
species_algo_str <- sprintf("%s_glm", trait.species)

# Geographically constrained modelling and merge the environmental data into trait.data
if (!is.null(trait.data)) {
  merged.result <- bccvl.trait.constraint.merge(trait.data, trait.data.params,
                                                environ.rasterstack,
                                                enviro.data.constraints,
                                                enviro.data.generateCHull)
  trait.data <- merged.result$data
  trait.data.params <- merged.result$params
  environ.constrained.rasterstack <- merged.result$raster
}

## MODEL

# Generate a formula for each trait
formulae_glm <- bccvl.trait.gen_formulae(trait.data.params)
for (formula in formulae) {
  
  # polr function for ordinal traits, multinom function for nominal traits, glm function for continuous traits
  na_action <- get0(EC.params$na_action)
  if (is.null(na_action)) {
    na_action = get("na.fail")
  }
  # Replace underscore and white-space with dash in trait-name for use in output filename
  traitname <- gsub("[_ ]", "-", trimws(formula$trait))
  if (formula$type == 'ordinal') {
    library(visreg)
    output_filename = paste0(traitname, "_polr_results.txt")
    glm.result <- polr(formula=formula(formula$formula),
                      data=trait.data,
                      weights=NULL,
                      na.action=na_action,
                      contrasts=NULL,
                      Hess=TRUE,
                      model=TRUE,
                      method="logistic")
    visreg(glm.result)
  } else if (formula$type == 'nominal') {
    output_filename = paste0(traitname, "_nom_results.txt")
    glm.result = multinom(formula=formula(formula$formula),
                          data=trait.data,
                          weights=NULL,
                          na.action=na_action,
                          contrasts=NULL,
                          summ=0,
                          model=TRUE)
  } else {
    output_filename = paste0(traitname, "_glm_results.txt")
    glm.result <- glm(formula=formula(formula$formula),
                     family=EC_FamilyFromString(EC.params$family),
                     data= trait.data,
                     weights=NULL,
                     na.action=na_action,
                     start=NULL,
                     etastart=NULL,
                     mustart=NULL,
                     offset=NULL,
                     model=TRUE,
                     method=EC.params$method,
                     x=FALSE,
                     y=FALSE,
                     contrasts=NULL)
    # regression and iagnostic plots
    bccvl.regPlots(glm.result, 
                   outerTitle = paste0(trait.species, ': GLM fit for ', formula$trait), 
                   fnamePrefix = paste0(traitname, '_', trait.species, '_glm_'))
    
    # Do projection only if there is no fixed factors i.e. all env variables of env dataset.
    env.names <- names(environ.rasterstack)
    if (all(unlist(lapply(formula$env, function(x) x %in% env.names)))) {
      # Do projection over constrained region
      if (is.null(environ.constrained.rasterstack)) {
        env.data <- as.data.frame(environ.rasterstack, xy=TRUE)
      }
      else {
        env.data <- as.data.frame(environ.constrained.rasterstack, xy=TRUE)
      }
      proj  <- predict.glm(glm.result, env.data)
      
      # combine the coordinates and projection value to generate a raster
      proj1 <- cbind(env.data[,c("x", "y")], proj)
      proj.raster <- rasterFromXYZ(proj1)
      
      # Save projection as geotif and png
      bccvl.saveModelProjection(proj.raster, projection.name, formula$trait, species_algo_str)
      
      # Do projection over unconstrained region only if all env variables are not categorical
      if (!is.null(environ.constrained.rasterstack) && all(unlist(lapply(environ.rasterstack@layers,
                                                                         function(lyr)
                                                                           return(!lyr@data@isfactor))))) {
        env.data = as.data.frame(environ.rasterstack, xy=TRUE)
        proj  = predict.glm(glm.result, env.data)
        
        # combine the coordinates and projection value to generate a raster
        proj1 = cbind(env.data[,c("x", "y")], proj)
        proj.raster = rasterFromXYZ(proj1)
        
        # Save projection as geotif and png
        bccvl.saveModelProjection(proj.raster, projection.name, formula$trait,s
                                  pecies_algo_str, filename_ext="unconstrained")
      }
    }
  }
  
  # Save the model to file
  EC_Save(glm.result, paste0(traitname, "_glm_model.object.RData"))
  
  ## Save the results as text to file for each trait
  s <- summary(glm.result) 
  bccvl.write.text(s, output_filename)
}


#===================================================================
## Trai Difference- GLM
##

## Runs a Generalized Linear Model to test how traits differ among species

# Load the library
library("MASS")
library("nnet")  

# Empty environmental dataset as it is not needed
environ.rasterstack = stack()
crs(environ.rasterstack) <- '+init=epsg:4326'

# Geographically constrained modelling; just need to constraint the trait-data
if (!is.null(trait.data)) {
  merged.result <- bccvl.trait.constraint.merge(trait.data, trait.data.params,
                                               environ.rasterstack,
                                               enviro.data.constraints,
                                               generateCHull=FALSE,
                                               generateGeoconstraint=FALSE)
  trait.data <- merged.result$data
  trait.data.params <- merged.result$params
}


## MODEL


# Generate the formulae to test differences among species for each trait separately
formulae_traitdiffglm <- bccvl.trait.gen_formulae(trait.data.params, trait_diff=TRUE)
# polr function for ordinal traits, multinom function for nominal traits, glm function for continuous traits
na_action = get0(EC.params$na_action)
if (is.null(na_action)) {
  na_action = get("na.fail")
}
for (formula in formulae) {
  traitname = gsub("[_ ]", "-", trimws(formula$trait))
  if (formula$type == 'ordinal') {
    output_filename = paste0(traitname, "_diffpolr_results.txt")
    glm.result <- polr(formula=formula(formula$formula),
                      data=trait.data,
                      weights=NULL,
                      na.action=na_action,
                      contrasts=NULL,
                      Hess=TRUE,
                      model=TRUE,
                      method="logistic")
  } else if (formula$type == 'nominal') {
    output_filename = paste0(traitname, "_diffnom_results.txt")
    glm.result <- multinom(formula=formula(formula$formula),
                          data=trait.data,
                          weights=NULL,
                          na.action=na_action,
                          contrasts=NULL,
                          summ=0,        
                          model=TRUE)
  } else {
    output_filename = paste0(traitname, "_diffglm_results.txt")
    glm.result <- glm(formula=formula(formula$formula),
                     family=family_from_string(EC.params$family),
                     data= trait.data,
                     weights=NULL,
                     na.action=na_action,
                     start=NULL,
                     etastart=NULL,
                     mustart=NULL,
                     offset=NULL,
                     model=TRUE,
                     method=EC.params$method,
                     x=FALSE,
                     y=FALSE,
                     contrasts=NULL)
  }
  
  ## Save the result to file
  # Save the model
  EC_Save(glm.result, paste0(traitname, "_diffglm_model.object.RData"))
  
  ## Save the results as text to file for each trait
  s <- summary(glm.result) 
  bccvl.write.text(s, output_filename)
}




#===================================================================
## GLMM
##

## DEFINITION

# Load the library
library(lme4)
library(ordinal)

# Generate a formula for each trait
# trait ~ fixed1 + fixed2 + (1|random1) + (1|random2)
formulae_glmm <- bccvl.trait.gen_formulae(trait.data.params, include_rf=TRUE)
for (formula in formulae) {
  # Run model - with clmm function for ordinal traits, glmer function for nominal traits, glmer function for continuous traits
  # Todo: not sure whether 'glmer' works for nominal trait data - need to further look into this
  na_action = get0(EC.params$na_action)
  if (is.null(na_action)) {
    na_action = get("na.fail")
  }
  
  traitname = gsub("[_ ]", "-", trimws(formula$trait))
  if (formula$type == 'ordinal') {
    output_filename = paste0(traitname, "_clmm_results.txt")
    glmm.result <- clmm(formula=formula(formula$formula),
                       data=trait.data,
                       weights=NULL,
                       na.action=na_action,
                       contrasts=NULL,
                       Hess=TRUE,
                       model=TRUE)
  } else if (formula$type == 'nominal') {
    output_filename = paste0(traitname, "_nom_results.txt")
    glmm.result <- glmer(formula=formula(formula$formula),
                        data=trait.data,
                        weights=NULL,
                        na.action=na_action,
                        contrasts=NULL,
                        summ=0,
                        model=TRUE)
  } else {
    output_filename = paste0(traitname, "_glmer_results.txt")
    glmm.result <- glmer(formula=formula(formula$formula),
                        family=EC_FamilyFromString(EC.params$family),
                        data= trait.data,
                        weights=NULL,
                        na.action=na_action,
                        start=NULL,
                        etastart=NULL,
                        mustart=NULL,
                        offset=NULL,
                        contrasts=NULL)
  }
  
  ## Save the result to file
  # Save the model
  EC_Save(glmm.result, paste0(traitname, "_glmm_model.object.RData"))
  
  ## Save the results as text to file for each trait
  s <- summary(glmm.result) 
  bccvl.write.text(s, output_filename)
}

