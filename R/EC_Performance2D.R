#' Calculate 2D measures of predictive performance for any
#' model that predicts a probability of presence (or success),
#' to be compared to binary observations of presence/absence (or
#' success/failure).
#'
#' @param obs 
#' @param pred 
#' @param species_algo_str 
#' @param make.plot 
#' @param kill.plot 
#'
#' @export
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @importFrom reshape2 melt
#' 
#' 


absmean <- function(x) abs(mean(x, na.rm=T))
absdiff <- function(x) abs(diff(x, na.rm=T))
pick1min <- function(x) { # w<-c(5,6,7,11,13)
  w <- which(x == min(x, na.rm=T));
  len <- length(w);
  if (len == 1) return(w) else {if (len %% 2 ==1) {return(median(w))} else {return(median(w[-1]))}}
}



EC_Performance2D <- function(obs, pred, species_algo_str, make.plot="EcoCommons", kill.plot=T) {
  library(gridExtra)
  library(pROC)

  # AIM: Calculate 2D measures of predictive performance for any
  # model that predicts a probability of presence (or success),
  # to be compared to binary observations of presence/absence (or
  # success/failure).
  #
  # AUTHOR: S.Low-Choy, Jan 2016, Griffith University
  # ACKNOWLEDGEMENTS: Shawn Laffan for useful discussions,
  # Chantal Huijbers, Sarah Richmond and Linda for testcase
  #
  # INPUTS
  # obs = vector of observations of presence/absence (binary)
  # pred = vector of predicted probabilities of presence
  #
  # OUTPUTS
  # Predicted presence/absence density plot and histogram
  # Sensitivity/Specificity plot
  # ROC plot
  # Plot with four different error rates across varying treshold probability values
  # 4 different loss functions (maximizing TPR+TNR, balancing all errors, or equalising error rates)
  # Table with best probability threshold value corresponding to minimum loss
  # Range of threshold probability values for which loss falls within 5% of minimum
  #

  # TESTING 23 July 2018 - SLC - problem with threshold = 0 when predictions only 0 or 1
  # species_algo_str <- "SRE"; make.plot <- "EC"; kill.plot <- F; obs<-gobs; pred<-gpred


  # TESTING 23 July 2018 - SLC - when more than one "best" value as the minimum is achieved twice!
  # Addition Number 2


  #############################################################
  #
  # ERROR CHECKING
  #
  # Check that observations and predictions match length
  nobs <- length(obs)
  if (nobs != length(pred)) stop("Ensure that vectors for observations and 
                                 predictions are of equal length!")

  # Check observations are binary, then recode as a factor if necessary
  tbl.obs <- table(obs)
  diversity.obs <- length(tbl.obs)
  if (diversity.obs==1) stop("All observations are the same!")
  if (diversity.obs!=2) stop("Ensure that the vector of observations only has two 
                             possible values.")

  # Check coding of presence and absence
  if (is.factor(obs)) {
    truth <- as.character(obs)
    temp <- as.numeric(truth)
    if (!is.na(temp))
      truth <- temp
  } else {
    truth <- obs
  }
  if (is.numeric(truth)) {
    if (any(truth==0)) {
      # presume 0 is coded to be absence
      # if presences are not 1s then recode them to be 1s
      if (any(truth[truth!=0]!=1)) {
        truth[truth!=0] <- 1
      }
    } else if (any(truth ==1)) {
      # if there are no zeros, then presume 1-2 coding
      if (any(truth ==2)) {
        truth <- truth-1
      }
    } else {
      stop("Can't figure out coding of absence-presence as it is not 0-1 or 1-2. 
           Suggest you recode obs as a factor, with levels 0 and 1")
    } # end if any truth == 0
  } else if (is.character(truth)) {
    # look for "p" for presence, "d" for detected or "o" for occupied
    the.pres <- grep("^[pP]", truth)
    if (any(the.pres)) {
      letter.pres <- "p"
    } else {
      the.pres <- grep("^[Dd]", truth)
      if (any(the.pres)) {
        letter.pres <- "d"
      } else {
        the.pres <- grep("^[oO]", truth)
        if (any(the.pres)) {
          letter.pres <- "o"
        } else {
          stop("Can't figure out coding of presences as they do not start with 
               the letter p, d or o. Suggest you recode presences using one of these options.")
        }
      }
    } # end if any the.pres
    truth[the.pres] <- 1

    # recode all non "p" or "d" to be absence (could be "a" for absence, or "n" for 
    # not present/seen/detected/occupied, or "u" for unseen etc, or "b" for background)
    tbl.obs <- table(truth)
    letter.abs <- names(tbl.obs)[-grep("1", names(tbl.obs))]
    truth[grep(paste("^",letter.abs,sep=""), truth)] <- 0

    if (any(the.pres)) {
      truth <- as.numeric(truth)
    }
  } else {
    stop("Ensure that the data type of observations is numeric, character or factor.")
  } # end if is.numeric(truth) or is.character(truth)

  # Check predictions are probabilities between 0 and 1
  if (any(pred < 0) | any(pred > 1)) stop("Predictions should be probabilities 
                                          between zero and one. (Check that predictions 
                                          are not on the log odds scale.)")

  #
  # MEASURES
  #

  # CREATE ERROR MATRICES
  list.tpv <- sort(unique(pred))

  tbl.pred <- table(pred)
  diversity.pred <- length(tbl.pred)
  if (diversity.pred==1) stop("All predictions are the same!")
  if (any(list.tpv==0)) { if (diversity.pred==2) { list.tpv <- 1; 
  } else { 
    list.tpv <- list.tpv[list.tpv!=0]; } }

  # tpv = threshold probability value: threshold that is used to transform the 
  # continuous probability of presence predicted by the model into a binary prediction:
  # probabilities < tpv = absences, probabilities > tpv = presences

  tp <- fp <- tn <- fn <- rep(NA, length(list.tpv))
  tpr <- fpr <- tnr <- fnr <- rep(NA, length(list.tpv))
  ppv <- npv <- fdr <- fors <- rep(NA, length(list.tpv))
  acc <- mcr <- tss <- bs <- rep(NA, length(list.tpv))
  tp.rand <- ets <- or <- csi <- rep(NA, length(list.tpv))
  Po <- Pe <- kappa <- roc <- auc <- rep(NA, length(list.tpv))

  # CALCULATE 1D MEASURES OF PREDICTIVE PERFORMANCE

  ## Elements of contigency table:

  # tp = True Positives (observed and predicted presences)
  # fp = False Positives (observed absences predicted as presences)
  # tn = True Negatives (observed and predicted absences)
  # fn = False Negatives (observed presences predicted as absences)

  for (ell in seq(along=list.tpv)) {  # ell <- 1

    th <- list.tpv[ell]

    tp[ell] <- length(which(pred>=th & truth==1))
    fp[ell] <- length(which(pred>=th & truth ==0))
    tn[ell] <- length(which(pred<th & truth ==0))
    fn[ell] <- length(which(pred<th & truth ==1))

    ## Evaluation statistics:

    # tpr = True Positive Rate (= Sensitivity) = proportion of observed presences that are correctly predicted.
    tpr[ell] <- tp[ell]/(tp[ell]+fn[ell])

    # fpr = False Positive Rate = proportion of observed absences that are incorrectly predicted.
    fpr[ell] <- fp[ell]/(fp[ell]+tn[ell])

    # tnr = True Negative Rate (= Specificity) = proportion of observed absences that are correctly predicted.
    tnr[ell] <- tn[ell]/(tn[ell]+fp[ell])

    # fnr = False Negative Rate = proportion of observed presences that are incorrectly predicted.
    fnr[ell] <- fn[ell]/(fn[ell]+tp[ell])

    # Note: tpr + fnr = 1, tnr + fpr = 1

    # ppv = Positive Predictive Value = for all predicted presences, how many were true observed presences?
    ppv[ell] <- tp[ell]/(tp[ell]+fp[ell])

    # npv = Negative Predictive Value = for all predicted absences, how many were true absences?
    npv[ell] <- tn[ell]/(tn[ell]+fn[ell])

    # fdr = False Discovery Rate = for all predicted presences, how many were false presences (observed absences)?
    fdr[ell] <- fp[ell]/(tp[ell]+fp[ell])

    # fors = False Omission Rate = for all predicted absences, how many were false absences (observed presences)?
    fors[ell] <- fn[ell]/(tn[ell]+fn[ell])

    # Note: ppv + fdr = 1, npv + fors = 1

    # acc = Accuracy = proportion of correctly predicted cases.
    acc[ell] <- (tp[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])

    # mcr = Misclassification Rate = proportion of incorrectly predicted cases.
    mcr[ell] <- (fp[ell]+fn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])

    # Note: acc + mcr = 1

    # tss = True Skill Statistic - describes how well the model separates presences from absences.
    tss[ell] <- (tpr[ell]-fpr[ell])

    # bs = Bias Score = frequency of predicted presences compared to the frequency of observed presences.
    bs[ell] <- (tp[ell]+fp[ell])/(tp[ell]+fn[ell])

    # csi = Critical Success Index = proportion of observed and predicted presences that are correct.
    csi[ell] <- tp[ell]/(tp[ell]+fp[ell]+fn[ell])

    # ets = Equitable Threat Score = proportion of observed and predicted presences that are correct, 
    # adjusted for true positives with random chance.
    tp.rand[ell] <- ((tp[ell]+fn[ell])*(tp[ell]+fp[ell]))/(tp[ell]+fp[ell]+tn[ell]+fn[ell])
    ets[ell] <- (tp[ell]-tp.rand[ell])/(tp[ell]+fp[ell]+fn[ell]-tp.rand[ell])

    # or = Odds Ratio = ratio of a correct prediction to an incorrect prediction.
    or[ell] <- (tp[ell]*tn[ell])/(fp[ell]*fn[ell])

    # kappa = Cohen's Kappa = Accuracy of the prediction relative to that of random chance.
    Po[ell] <- acc[ell] # observed accuracy
    Pe[ell] <- ((((tp[ell]+fp[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))*((tp[ell]+fn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))) + (((fn[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell]))*((fp[ell]+tn[ell])/(tp[ell]+fp[ell]+tn[ell]+fn[ell])))) # expected accuracy by random chance
    kappa[ell] <- ((Po[ell]-Pe[ell])/(1-Pe[ell]))
  }

  # auc = Area Under the (ROC) Curve
  roc <- pROC::roc(truth, pred)
  auc <- pROC::auc(roc)

  # Compile the information into dataframes
  temp <- data.frame(list(tpv=list.tpv, tpr=tpr, fpr=fpr,  tnr=tnr,
                          fnr=fnr, ppv=ppv, fdr=fdr, npv=npv, fors=fors))
  auc.d <- data.frame(auc)
  auc.d <- round(auc.d, digits = 2)

  # CALCULATE 2D MEASURES OF PREDICTIVE PERFORMANCE

  # Calculate losses as cost functions across errors
  list.errors <- c("fpr","fnr","fors","fdr")

  # Loss functions:
  # L.diag = minimizing the sum of the diagnostic errors (FPR, FNR)
  #        = maximizing the sum of the True Positive Rate (Sensitivity) and the True Negative Rate (Specificity)
  # L.pred = minimizing the sum of the predictive errors (FDR, FORS)
  #        = maximizing the sum of the Positive Predictive Value and the Negative Predictive Value
  # L.all = balancing all error rates: FPR, FNR, FDR, FORS
  # L.eq.diag = equalising diagnostic errors: FPR = FNR = cross-over of TPR and TNR

  temp$L.diag <- apply(temp[, c("fpr","fnr")],1, absmean)
  temp$L.pred <- apply(temp[, c("fors","fdr")],1, absmean)
  temp$L.all <- apply(temp[, list.errors],1, absmean)
  temp$L.eq.diag <- apply(temp[, c("tpr", "tnr")], 1, absdiff)

  # Addition Number 2, 23 July 2018
  # Check if there is more than one minimum, then pick the middle
  best <- list(
    diag = temp$tpv[pick1min(temp$L.diag) ],
    pred = temp$tpv[pick1min(temp$L.pred == min(temp$L.pred, na.rm=T))],
    all = temp$tpv[pick1min(temp$L.all == min(temp$L.all, na.rm=T))],
    eq.diag = temp$tpv[pick1min(temp$L.eq.diag==min(temp$L.eq.diag, na.rm=T))]
  )

  # End Addition Number 2, 23 July 2018

  # Calculate the range of threshold probability values for which each of 
  # the losses fall within 5% of the best value
  rangeperf <- matrix(NA, nrow=4, ncol=2, dimnames=list(names(best), c("lower","upper")))
  for (v in names(best)) { # v<-"eq.diag"
    the.v <- paste("L.", v, sep="")
    min.v <- min(temp[,the.v], na.rm=T)
    d.minv <- (abs(temp[,the.v] - min.v) / min.v)
    the.range <- temp$tpv[ d.minv < 0.05 ]
    rangeperf[dimnames(rangeperf)[[1]]==v, 1:2] <- c(min(the.range, na.rm=T), max(the.range, na.rm=T))
  }

  rangeperf <- as.data.frame(rangeperf)
  rangeperf$type.of.loss <- names(best)
  rangeperf$best <- unlist(best)

  loss.table <- subset(rangeperf, select = c("lower", "upper", "best"))
  row.names(loss.table) = c("Maximize TPR+TNR", "Maximize PPV+NPV", "Balance all errors", "TPR = TNR")
  loss.table <- as.data.frame(loss.table)

  # Rescale
  temp$L.eq.diag <- temp$L.eq.diag/max(temp$L.eq.diag, na.rm=T)


  #########################################################################
  ### 3: Functions to create SDM outputs
  #########################################################################

  if (make.plot!="") {
    # reshape the data so that it is in long rather than wide format
    # (= each row represents one item, labels are specified by 'measure' column; used by ggplot2)
    errs <- reshape2::melt(temp, id.var="tpv", measure.var=c("tpr", "tnr", "fpr", "fnr", "fdr",
                                                   "fors", "L.diag", "L.pred", "L.all", "L.eq.diag"))
    names(errs)[2] <- c("measure")

    # Create Presence/absence density plot across threshold probability values
    temp2 <- data.frame(list(pred=pred, obs=obs))
    png(file=file.path(EC.env$outputdir, sprintf("%s-presence-absence-plot_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g1 <- ggplot(temp2, aes(x=pred, fill=factor(obs))) +
      geom_density(stat="density", alpha=0.5) +
      labs(title="Presence/absence density plot \nacross predicted probability of presence",
           x="\nPredicted probability of presence", y="Density\n") +
      scale_fill_manual(values=c("#EE3B3B", "#6495ED"), labels=c(" Absences      ", " Presences")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(), legend.key.size=unit(1.5, "lines")) +
      guides(fill=guide_legend(nrow=1, title=NULL))
    print(g1)
    dev.off()

    # Create Presence/absence histogram across threshold probability values
    png(file=file.path(EC.env$outputdir, sprintf("%s-presence-absence-hist_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g2 <- ggplot(temp2, aes(x=pred, fill=factor(obs)))  +
      geom_histogram(position="dodge", alpha = 0.5) +
      labs(title="Presence/absence histogram \nacross predicted probability of presence",
           x="\nPredicted probability of presence", y="Count\n") +
      scale_fill_manual(values=c("#EE3B3B", "#6495ED"), labels=c(" Absences    ", " Presences")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(), legend.key.size=unit(1.5, "lines")) +
      guides(fill=guide_legend(nrow=1, title=NULL))
    print(g2)
    dev.off()

    # Create TPR-TNR plot
    png(file=file.path(EC.env$outputdir, sprintf("%s-TPR-TNR_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g3 <- ggplot(errs[errs$measure %in% c("tpr", "tnr"), ],
                 aes(x=tpv, y=value, colour=measure)) +
      geom_line(size=1.2) +
      ylim(0,1) +
      labs(title="Sensitivity-Specificity plot\n", x="\nThreshold probability value",
           y="TPR/TNR value\n") +
      scale_colour_manual(values=c("#3CAB34", "#049CE3"), labels=c("True Positive Rate (=Sensitivity)",
                                                                   "True Negative Rate (=Specificity)")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(),
            legend.key.size=unit(2.5, "lines")) +
      guides(colour=guide_legend(nrow=2, title=NULL))
    print(g3)
    dev.off()

    # Create Error rates plot: shows the values of four different error rates across the range of threshold probability values
    png(file=file.path(EC.env$outputdir, sprintf("%s-error-rates_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g4 <- ggplot(errs[errs$measure %in% c("fpr", "fnr", "fdr", "fors"), ],
                 aes(x=tpv, y=value, colour=measure, linetype=measure)) +
      geom_line(size=1.2) +
      ylim(0,1) +
      labs(title="Error rates plot\n", x="\nThreshold probability value", y="Error rate value\n") +
      scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed"),
                            labels=c("False Positive Rate  ", "False Negative Rate   ",
                                     "False Discovery Rate", "False Omission Rate")) +
      scale_colour_manual(values=c("#FAB334", "#D55E00", "#FAB334", "#D55E00"),
                          labels=c("False Positive Rate  ", "False Negative Rate   ",
                                   "False Discovery Rate", "False Omission Rate")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(),
            legend.key.size=unit(2.5, "lines")) +
      guides(colour=guide_legend(nrow=2, title=NULL), linetype=guide_legend(nrow=2, title=NULL))
    print(g4)
    dev.off()

    # Create ROC plot
    png(file=file.path(EC.env$outputdir, sprintf("%s-ROC_%s.png", make.plot,
                                                 species_algo_str)), width=480, height=480)
    xmax1 = min(round((max(temp$fpr) + 0.2)/0.1)*0.1, 1)
    xpos = max(xmax1/2, 0.1)
    g5 <- ggplot(temp, aes(x=fpr, y=tpr)) +
      geom_line(size=1.2) +
      ylim(0,1) +
      xlim(0, xmax1) +
      geom_abline(intercept=0, slope=1, colour="grey") +
      labs(x="\nFalse Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)\n") +
      ggtitle(paste("ROC plot")) +
      annotate(geom = "text", x = xpos, y = 0.1, label = paste("AUC = ", auc.d$auc), size = 6) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)))
    print(g5)
    dev.off()

    # Create evaluation stats table with values for optimum tpv value.
    all.stats <- data.frame(list(tpv=list.tpv, tpr=tpr, tnr=tnr, fpr=fpr,
                                 fnr=fnr, fdr=fdr, fors=fors, ppv=ppv, npv=npv,
                                 kappa=kappa, tss=tss, bs=bs, csi=csi, ets=ets,
                                 or=or, acc=acc, mcr=mcr))
    all.stats <- round(all.stats, digits = 3)

    # Addition Number 2, 23 July 2018
    selected_rows = c(pick1min(temp$L.diag),
                      pick1min(temp$L.pred == min(temp$L.pred, na.rm=T)),
                      pick1min(temp$L.all == min(temp$L.all, na.rm=T)),
                      pick1min(temp$L.eq.diag == min(temp$L.eq.diag, na.rm=T)))
    stats.table <- all.stats[selected_rows, ] # select row with all stats for maximum value of loss methods
    rownames(stats.table) <- c("max TPR + TNR", "max PPV + NPV", "balance all errors", "TPR=TNR")
    # End Addition Number 2, 23 July 2018

    # stats.table <- rbind(max.TPR.TNR, TPR.eq.TNR, max.Kappa)
    names(stats.table) <- c("Optimum threshold value:", "True Positive Rate (TPR)",
                            "True Negative Rate (TNR)", "False Positive Rate (FPR)",
                            "False Negative Rate (FNR)", "False Discovery Rate (FDR)",
                            "False Omission Rate (FOR)", "Positive Predictive Value (PPV)",
                            "Negative Predictive Value (NPV)", "Cohen's Kappa",
                            "True Skill Statistic (TSS)", "Bias Score (BS)",
                            "Critical Success Index (CSI)", "Equitable Threat Score (ETS)",
                            "Odds-Ratio (OR)","Accuracy", "Misclassification Rate")
    eval.stats <- t(stats.table) # transpose table

    # Create Loss function plot: shows the values of different loss functions
    # across the range of threshold probability values
    png(file=file.path(EC.env$outputdir, sprintf("%s-loss-functions_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g6 <- ggplot(errs[errs$measure %in% rev(c("L.diag", "L.pred", "L.all", "L.eq.diag")), ],
                 aes(x=tpv, y=value, colour=measure)) +
      geom_line(size=1.2) +
      ylim(0,1) +
      labs(title="Loss function plot\n", x="\nThreshold probability value",
           y="Loss function value\n") +
      scale_colour_manual(values=c("#48D1CC", "#9F79EE", "#EE9572", "#FF3E96"),
                          labels=c("Maximize TPR + TNR   ", "Maximize PPV + NPV   ",
                                   "Balance all errors", "TPR = TNR")) +
      theme(axis.text = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(),
            legend.key.size=unit(2.5, "lines")) +
      guides(colour=guide_legend(nrow=2, title = NULL))
    print(g6)
    dev.off()

    # Create Loss functions-intervals plot within 5% of the best value
    rangeperf$type.of.loss <- factor(rangeperf$type.of.loss,
                                     levels=(c("diag", "pred", "all", "eq.diag")))
    png(file=file.path(EC.env$outputdir, sprintf("%s-loss-intervals_%s.png",
                                                 make.plot, species_algo_str)), width=480, height=480)
    g7 <- ggplot(rangeperf, aes(x=type.of.loss, y=best, ymin=lower, ymax=upper,
                                colour=type.of.loss)) +
      geom_pointrange(size=1.2) +
      geom_line(size=1.2) +
      coord_flip() +
      ylim(0,1) +
      scale_x_discrete(limits=c("diag","pred", "all", "eq.diag")) +
      scale_colour_manual(values=c("#48D1CC", "#9F79EE", "#EE9572", "#FF3E96"),
                          labels=c("Maximize TPR + TNR   ", "Maximize PPV + NPV   ",
                                   "Balance all errors   ", "TPR = TNR")) +
      labs(title="Range of threshold probability value \nwithin 5% of minimum per loss\n",
           x="Type of loss function\n", y="\nThreshold probability value") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x = element_text(family="Arial", size=rel(1.5)),
            axis.title = element_text(family="Arial", size=rel(1.5)),
            plot.title = element_text(family="Arial", size=rel(2)),
            legend.text = element_text(family="Arial", size=rel(1.5)),
            legend.position="top", legend.key=element_blank(),
            legend.key.size=unit(2.5, "lines")) +
      guides(colour=guide_legend(nrow=2, title = NULL))
    print(g7)
    dev.off()
  }

  return(list(performance=temp, stats=eval.stats, loss.summary=loss.table))
}
