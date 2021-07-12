EC_VIPplot <- function(fittedmodel=NULL,
                       method=c("glm","cta","gam","ann", "rf", "gbm", "mars", "maxent"),
                       cor.method=c("pearson","spearman"),
                       pdf=TRUE, biom_vi=FALSE,output.table=FALSE, data1, this.dir, filename)
{
  library("ggplot2")
  library("reshape2")
  library("mgcv")
  library("rpart")
  library("caret")
  library("ggdendro")

  # README notes:
  # (1) fittedmodel: the fitted model object obtained running the biomod2 function'BIOMOD_Modeling'.
  #     method is one of "glm","rpart","gam","ann", "rf", "gbm", "mars", "maxent.p"
  #     cor.method is one of "pearson","spearman"
  # (2) pdf=FALSE: the default setting - no pdf file generated;  pdf=TRUE: the VIP generated is
  #     saved as a pdf file in the working directory.
  # (3) biom_vi=FALSE: a function/algorithm other than the biomod2 inbuilt function 'variables_importance'
  #     will be applied for evaluating/ranking the variable importance.
  # (4) output.table=FALSE: a .csv file which contains a table to display the glm model parameter
  #     estimates and the 95% confidence bounds for both raw data and scaled data will be generated
  #     if output.table=TRUE.
  # (5) cor.method=c("pearson","spearman"): the default "pearson" method measures only the linear
  #     association among variables; the "spearman" method is a rank-based algorithm which is more
  #     robust measure for association among variables, e.g., the non-linear association will also be detected.
  # (6) Note that 'data1' is a dataframe with the response variable in the first column following
  #     by the predictor variables / enviornmental variables in the other columns;
  #     'data1' is needed for generating the VIP by the biomod2 inbuilt function
  #     'variables_importance(fittedmodel, data1)' and for calculate the AIC scores.
  #     Warning: Except for the glm algorithm, the 'data1' should include all predictor variables to make
  #     the variable importance ranking outcomes meaningful.
  # (7) this.dir specifies the route to access the biomod2 model (but not including the model name)
  # (8) filename to be saved without the file extension.
  #

  data1$y.data[is.na(data1$y.data)] <- 0

  # extract the root of filenames used by biomod to save model results
  filenames <- dir(this.dir)
  #loc <- regexpr("_RUN[[:digit:]]", filenames[1])
  #fileroot <- substr(filenames[1], 1, loc-1)

  # select the full model generated
  filekeep <-  paste(this.dir, "/", filenames[1], sep="")

  if (!is.na(match("glm",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    fitteddata = fittedmodel$data  # these are the data used by biomod2 for model fitting
    nd = dim(fitteddata)[2]
    sub.data = fitteddata[,2:nd, drop=FALSE]
    # Can only scale numeric data
    cat.data = Filter(is.factor, sub.data)
    sub.data = Filter(is.numeric, sub.data)
    RespV = fitteddata[,1, drop=FALSE]
    rescaled.data <- scale(sub.data)

    # attributes(rescaled.data)
    all.data <- cbind(RespV, as.data.frame(rescaled.data), cat.data)

    # head(all.data); dim(all.data)
    rescaled.glm <- glm(fittedmodel$formula, data=all.data, family=binomial)

    scaled = as.data.frame(cbind(coef(rescaled.glm), confint(rescaled.glm)))
    scaled = scaled[-1,]
    raw = as.data.frame(cbind(coef(fittedmodel), confint(fittedmodel)))
    raw = raw[-1,]

    #  variable importance plot in terms of relative effect size
    #  the relative effect size is defined as the coefficient estimated based on scaled data
    nx = length(coef(fittedmodel)[-1])
    df1 = as.data.frame(cbind(1:nx,round(scaled,3)))

    names(df1) = c("xvar", "meanest", "lower", "upper")
    df1$xvar = factor(df1$xvar, labels = rownames(df1))

    p1 <- ggplot(df1, aes(x=xvar, y=meanest)) + geom_hline(yintercept = 0) + labs(x=" ")

    ps = p1 + geom_errorbar(aes(ymin=lower,ymax=upper),lwd=0.8,width=0.25) +
      labs(y="relative effect size") +  labs(title="         scaled data") + coord_flip()

    df2 = as.data.frame(cbind(1:nx,round(raw,3)))
    names(df2) = c("xvar", "meanest", "lower", "upper")
    df2$xvar = factor(df2$xvar, labels = rownames(df2))

    if (output.table)
    {
      df1t = df1; df2t = df2
      names(df1t) = c("x.var", "coeff.est.scaled", "lower", "upper")
      names(df2t) = c("x.var", "coeff.est.raw", "lower", "upper")
      dfout = cbind(df2t,df1t)
      write.csv(dfout,file=paste(filekeep,"paraest_out.csv",sep="_"),row.names=FALSE)
    }

    #  the heatmap in terms of correlation among numerical predictor variables
    rescdata = Filter(is.numeric, rescaled.glm$data[,-1, drop=FALSE])

    if("spearman" %in% cor.method) {
      xx = cor(rescdata, method="spearman")
    } else if("pearson" %in% cor.method)  {
      xx = cor(rescdata)
    }

    lower_tri <- xx
    lower_tri[upper.tri(lower_tri)] <- NA

    xx.ml <- melt(lower_tri,na.rm=TRUE)  #the argument 'na.rm=TRUE' seems not working)

    corx = xx.ml[,3]
    rm = which(is.na(corx)==TRUE)
    xx.ml = xx.ml[-rm,]

    pheat <- ggplot(xx.ml, aes(X1, X2)) + geom_tile(aes(fill = value), colour="black") +
      scale_fill_gradient2(low = "green4", high = "violetred", mid="white",
                           midpoint=0, limit=c(-1,1)) + labs(y=" ") + theme_minimal() +
      scale_x_discrete(limits=rownames(xx)) + scale_y_discrete(limits=colnames(xx)) + coord_fixed() +
      theme(axis.title.x=element_blank(),legend.position = "bottom", axis.text.x=element_text(angle=-90)) +
      guides(fill=guide_legend(title="correlation"))

    # Save as variable correlation plot.
    filename1 = sub("vip_plot", "variable_correlations", filename)
    EC_SavePDF(pheat, ncol=1, nrow=1, filename=filename1, aspdf=pdf)

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
    pa <- ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") +
      labs(y="AIC score for information loss") + labs(title="AIC approach")

    ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()

    # the variable importance plot using the inbuilt biomod2 function 'variables_importance'
    vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
    nx = length(vi_biomod)
    dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
    pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")

    ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

    # Save as variable relative contribution plot.
    filename1 = sub("vip_plot", "variable_relative_contribution", filename)
    if (biom_vi) {
      EC_SavePDF(ps, ppv, ncol=2, nrow=1, filename=filename1, aspdf=pdf)
    }
    else {
      EC_SavePDF(ps, ppa, ncol=2, nrow=1, filename=filename1, aspdf=pdf)
    }
  }

  if (!is.na(match("cta",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    # variable importance is part of the model fitting outcomes with 'rpart' algorithm and
    # this information can be used for generating the variable importance plot

    varimp0 = fittedmodel$variable.importance
    nx = length(varimp0)
    df0 = as.data.frame(cbind(1:nx,varimp0))
    df0$V1 = factor(df0$V1, labels = rev(names(varimp0)))

    p <- ggplot(df0, aes(x=V1, y=rev(varimp0))) + labs(x=" ") +
      labs(y="variable importance score") + labs(title="part of the 'rpart' model output")

    pp0 = p + geom_col(alpha=0.6,col="blue") + coord_flip()

    ddata <- dendro_data(fittedmodel)
    ppt = ggplot() +
      geom_segment(data = ddata$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_text(data = ddata$labels, aes(x = x, y = y, label = label), size = 3, vjust = 0) +
      geom_text(data = ddata$leaf_labels, aes(x = x, y = y, label = label), size = 3, vjust = 1) +
      theme_dendro()

    # variable importance plot using the inbuilt biomod2 function 'variables_importance'

    nd = dim(data1)[2]
    subdata1 = data1[,2:nd, drop=FALSE]
    vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
    nx = length(vi_biomod)
    dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
    pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
    ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

    if (biom_vi)
    {
      EC_SavePDF(ppt, ppv, ncol=2, nrow=1, filename=filename, aspdf=pdf)
    }
    else
    {
      EC_SavePDF(ppt, pp0, ncol=2, nrow=1, filename=filename, aspdf=pdf)
    }
  }


  if (!is.na(match("gam",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    if (biom_vi)
    {
      # variable importance plot using the inbuilt biomod2 function 'variables_importance'
      nd = dim(data1)[2]
      RespV1 = data1[,1]; subdata1 = data1[,2:nd, drop=FALSE]

      vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
      nx = length(vi_biomod)
      dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
      dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
      pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
        labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
      ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

      EC_SavePDF(ppV, ncol=1, nrow=1, filename=filename, aspdf=pdf)
    }  # end of 'if(biom_vi=TRUE)'
    else
    {
      # variable importance plot following the AIC approach
      nd = dim(data1)[2]
      RespV1 = data1[,1]
      subdata1 = data1[,2:nd, drop=FALSE]

      # gam function cannot take categorical data, so exclude categorical data.
      subdata1 = Filter(is.numeric, subdata1)

      xname = names(subdata1)
      sname = paste("s(", xname, ")",sep="")

      gamformu.all <- as.formula(paste("RespV1 ~ 1 +", paste(sname, collapse= "+")))
      gam.all = gam(formula = gamformu.all, family = binomial, data = subdata1)

      Xaic = NULL
      nd = dim(subdata1)[2]
      for (i in 1:nd)
      {
        subdf = subdata1[, -i, drop=FALSE]
        xname1 = names(subdf)
        sname1 = paste("s(", xname1, ")",sep="")
        gamformu1 <- as.formula(paste("RespV1 ~ 1 +", paste(sname1, collapse= "+")))
        gam.one = gam(formula = gamformu1, family = binomial, data = subdf)
        Xaic = c(Xaic,AIC(gam.one))
      }

      relaAIC = round(Xaic - AIC(gam.all),2)
      nx = length(relaAIC)
      dfa = as.data.frame(cbind(1:nx,relaAIC))
      dfa$V1 = factor(dfa$V1, labels = rev(names(subdata1)))
      pa <- ggplot(dfa, aes(x=V1, y=rev(relaAIC))) + labs(x="predictor variables") +
        labs(y="AIC score for information loss") + labs(title="AIC approach")

      ppa = pa + geom_col(alpha=0.6,col="blue") + coord_flip()

      EC_SavePDF(ppa, ncol=1, nrow=1, filename=filename, aspdf=pdf)

    }
  }

  if (!is.na(match("ann",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    # variable importance plot using the inbuilt biomod2 function 'variables_importance'
    nd = dim(data1)[2]
    RespV1 = data1[,1]; subdata1 = data1[,2:nd]

    vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
    nx = length(vi_biomod)
    dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
    pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
    ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

    EC_SavePDF(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  }

  if (!is.na(match("mars",method)))
  {
    # Note that the model class generated from MARS algorithm is not supported by
    #  the inbuilt biomod2 function 'variables_importance', neither the AIC approach is applicable.
    # However, the function 'varImp' in package 'caret' accept the MARS model object for estimating
    #  the variable importance. GCV = generalized cross validation
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    # variable importance plot using the inbuilt function 'varImp' from package 'caret'
    nd = dim(data1)[2]
    RespV1 = data1[,1]; subdata1 = data1[,2:nd]

    var_imp = varImp(fittedmodel)
    nx = length(var_imp[,1])
    dfvi = as.data.frame(cbind(1:nx,var_imp[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(var_imp)))
    pv <- ggplot(dfvi, aes(x=V1, y=rev(var_imp[,1]))) + labs(x="predictor variables") +
      labs(y="relative reduction in GCV") + labs(title="function 'varImp' in package 'caret'")
    ppv = pv + geom_col(alpha=0.6,col="red") + coord_flip()

    EC_SavePDF(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  }


  if (!is.na(match("gbm",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    # variable importance plot using the inbuilt biomod2 function 'variables_importance'
    nd = dim(data1)[2]
    RespV1 = data1[,1]; subdata1 = data1[,2:nd]

    vi_biomod = variables_importance(fittedmodel,data=subdata1)$mat
    nx = length(vi_biomod)
    dfvi = as.data.frame(cbind(1:nx,vi_biomod[,1]))
    dfvi$V1 = factor(dfvi$V1, labels = rev(rownames(vi_biomod)))
    pv <- ggplot(dfvi, aes(x=V1, y=rev(vi_biomod[,1]))) + labs(x="predictor variables") +
      labs(y="variable importance score") + labs(title="biomod2 function 'variables_importance'")
    ppv = pv + geom_col(alpha=0.6,col="green4") + coord_flip()

    EC_SavePDF(ppv, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  }


  if (!is.na(match("rf",method)))
  {
    working <- load(filekeep)
    fittedmodel <- get_formal_model(eval(parse(text=working)))

    # Random forests (rf) provide an improvement over bagged trees by way of a small tweak
    #  that decorrelates the trees.
    # Note that the variable importance plot using the inbuilt biomod2 function 'variables_importance'
    #  does not seem working with Random forests algorithm.  On the other hand, however, the fitted
    #  rf model object contains the variable importance information which is measured by the mean derease
    #  in Gini index (expressed relative to the maximum).
    # While RSS is used for measuring the regression tree model performance, the Gini index is used for
    #  measuring the classification tree model performance and Gini index is a measure of total variance
    #  across the K classes.

    nd = dim(data1)[2]
    RespV1 = data1[,1]; subdata1 = data1[,2:nd]

    out.rf = fittedmodel$importance
    rfImp = out.rf[,1]
    nx = length(rfImp)

    dfrf = as.data.frame(cbind(1:nx,rfImp))
    dfrf$V1 = factor(dfrf$V1, labels = rev(names(rfImp)))
    prf <- ggplot(dfrf, aes(x=V1, y=rev(rfImp))) + labs(x="predictor variables") +
      labs(y="mean decrease in Gini index") + labs(title="part of rf model fitting outputs")
    pprf = prf + geom_col(alpha=0.6,col="blue") + coord_flip()

    EC_SavePDF(pprf, ncol=1, nrow=1, filename=filename, aspdf=pdf)
  }


  if (!is.na(match("maxent",method)))
  {
    # variable importance is part of the model fitting outcomes with 'MAXENT.Phillips' algorithm
    if (regexpr("_outputs", filekeep) < 0)
    {
      working <- paste(filekeep,"_outputs/maxentResults.csv",sep="")
    }
    else
    {
      working <- paste(filekeep,"/maxentResults.csv",sep="")
    }

    df.P = read.csv(working)

    the.data <- data1[,-1]
    the.data <- Filter(is.numeric, the.data)
    nx = dim(the.data)[2]    #decide the number of the predictor variables
    yp = as.numeric(df.P[,8:(7+nx)])

    dfp = as.data.frame(cbind(1:nx,yp))
    dfp$V1 = factor(dfp$V1, labels = rev(names(the.data)))
    p <- ggplot(dfp, aes(x=V1, y=rev(yp))) + labs(x="predictor variables") +
      labs(y="variable relative contribution (%)") + labs(title="maxent algorithm")
    pp = p + geom_col(alpha=0.6,col="red") + coord_flip()

    # Save as variable relative contribution plot.
    filename1 = sub("vip_plot", "variable_relative_contribution", filename)
    EC_SavePDF(pp, ncol=1, nrow=1, filename=filename1, aspdf=pdf)

    #  the heatmap in terms of correlation among predictor variables
    if(cor.method=="pearson")  xx = cor(the.data)
    if(cor.method=="spearman") xx = cor(the.data,method="spearman")

    get_lower_tri<-function(cormat)
    {
      cormat[upper.tri(cormat)] <- NA
      return(cormat)
    }

    lower_tri = get_lower_tri(xx)

    xx.ml <- melt(lower_tri,na.rm=TRUE)  #the argument 'na.rm=TRUE' seems not working)

    corx = xx.ml[,3]
    rm = which(is.na(corx)==TRUE)
    xx.ml = xx.ml[-rm,]

    pheat <- ggplot(xx.ml, aes(X1, X2)) + geom_tile(aes(fill = value), colour="black") +
      scale_fill_gradient2(low = "green4", high = "violetred", mid="white", midpoint=0, limit=c(-1,1)) +
      scale_x_discrete(limits=rownames(xx)) + scale_y_discrete(limits=colnames(xx)) + coord_fixed() +
      labs(y=" ") + theme_minimal() +
      theme(axis.title.x=element_blank(),legend.position = "bottom", axis.text.x=element_text(angle=-90)) +
      guides(fill=guide_legend(title="correlation"))
    # Save as variable correlation plot.
    filename1 = sub("vip_plot", "variable_correlations", filename)
    EC_SavePDF(pheat, ncol=1, nrow=1, filename=filename1, aspdf=pdf)
  }
}
