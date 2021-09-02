#' Summary statistics to evaluate model performance
#' Return equal sensitivity specificity threshold
#'
#' @param occur 
#' @param absen 
#' @param model 
#' @param predictor 
#' @param species 
#'
#' @export EC_SummaryStats
#' @importFrom dismo evaluate

EC_SummaryStats <- function(occur, absen, model, predictor, species) {
  writeLines('Compute evaluation statistics ...')
  evaluation<-dismo::evaluate(occur, absen, model, predictor)  # statistics for different threshold methods
  
  png(file=file.path(EC.env$outputdir, paste0("density_", species, "_maxent.png")),
      width=800, height=800)
  density(evaluation)
  dev.off()
  
  # Retrieve threshold values for each method
  threshold_all <- threshold(evaluation)
  threshold_all <- subset(threshold_all[,1:5]) # leave out sensitivity threshold
  threshold_kappa <- threshold(evaluation, stat='kappa')
  threshold_spec_sens <- threshold(evaluation, stat='spec_sens')
  threshold_no_omission <- threshold(evaluation, stat='no_omission')
  threshold_prevalence <- threshold(evaluation, stat='prevalence')
  threshold_equal_sens_spec <- threshold(evaluation, stat='equal_sens_spec')
  
  # Run the evaluate function for each threshold method
  evaluation_kappa <- dismo::evaluate(occur, absen, model, predictor, threshold_kappa)
  evaluation_spec_sens <- dismo::evaluate(occur, absen, model, predictor, threshold_spec_sens)
  evaluation_no_omission <- dismo::evaluate(occur, absen, model, predictor, threshold_no_omission)
  evaluation_prevalence <- dismo::evaluate(occur, absen, model, predictor, threshold_prevalence)
  evaluation_equal_sens_spec <- dismo::evaluate(occur, absen, model, predictor, threshold_equal_sens_spec)
  
  # Create a dataframe of evaluation statistics using the threshold where Kappa is highest
  kappa_prevalence = t(as.data.frame(as.double(slot(evaluation_kappa,"prevalence"))))
  kappa_OR = t(as.data.frame(as.double(slot(evaluation_kappa,"OR"))))
  kappa_MCR = t(as.data.frame(as.double(slot(evaluation_kappa,"MCR"))))
  kappa_NPP = t(as.data.frame(as.double(slot(evaluation_kappa,"NPP"))))
  kappa_PPP = t(as.data.frame(as.double(slot(evaluation_kappa,"PPP"))))
  kappa_FNR = t(as.data.frame(as.double(slot(evaluation_kappa,"FNR"))))
  kappa_FPR = t(as.data.frame(as.double(slot(evaluation_kappa,"FPR"))))
  kappa_TNR = t(as.data.frame(as.double(slot(evaluation_kappa,"TNR"))))
  kappa_TPR = t(as.data.frame(as.double(slot(evaluation_kappa,"TPR"))))
  kappa_CCR = t(as.data.frame(as.double(slot(evaluation_kappa,"CCR"))))
  kappa_ODP = t(as.data.frame(as.double(slot(evaluation_kappa,"ODP"))))
  
  kappa_bind <- rbind(kappa_CCR, kappa_FNR, kappa_FPR, kappa_MCR, kappa_NPP,
                      kappa_ODP, kappa_OR, kappa_PPP, kappa_prevalence, kappa_TNR, kappa_TPR)
  
  
  # Create a dataframe of evaluation statistics using the threshold at which 
  # the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
  spec_sens_prevalence = t(as.data.frame(as.double(slot(evaluation_spec_sens,"prevalence"))))
  spec_sens_OR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"OR"))))
  spec_sens_MCR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"MCR"))))
  spec_sens_NPP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"NPP"))))
  spec_sens_PPP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"PPP"))))
  spec_sens_FNR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"FNR"))))
  spec_sens_FPR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"FPR"))))
  spec_sens_TNR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"TNR"))))
  spec_sens_TPR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"TPR"))))
  spec_sens_CCR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"CCR"))))
  spec_sens_ODP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"ODP"))))
  
  spec_sens_bind <- rbind(spec_sens_CCR, spec_sens_FNR, spec_sens_FPR, spec_sens_MCR,
                          spec_sens_NPP, spec_sens_ODP, spec_sens_OR, spec_sens_PPP,
                          spec_sens_prevalence, spec_sens_TNR, spec_sens_TPR)
  
  # Create a dataframe of evaluation statistics using the highest threshold at which there is no omission
  no_omission_prevalence = t(as.data.frame(as.double(slot(evaluation_no_omission,"prevalence"))))
  no_omission_OR = t(as.data.frame(as.double(slot(evaluation_no_omission,"OR"))))
  no_omission_MCR = t(as.data.frame(as.double(slot(evaluation_no_omission,"MCR"))))
  no_omission_NPP = t(as.data.frame(as.double(slot(evaluation_no_omission,"NPP"))))
  no_omission_PPP = t(as.data.frame(as.double(slot(evaluation_no_omission,"PPP"))))
  no_omission_FNR = t(as.data.frame(as.double(slot(evaluation_no_omission,"FNR"))))
  no_omission_FPR = t(as.data.frame(as.double(slot(evaluation_no_omission,"FPR"))))
  no_omission_TNR = t(as.data.frame(as.double(slot(evaluation_no_omission,"TNR"))))
  no_omission_TPR = t(as.data.frame(as.double(slot(evaluation_no_omission,"TPR"))))
  no_omission_CCR = t(as.data.frame(as.double(slot(evaluation_no_omission,"CCR"))))
  no_omission_ODP = t(as.data.frame(as.double(slot(evaluation_no_omission,"ODP"))))
  
  no_omission_bind <- rbind(no_omission_CCR, no_omission_FNR, no_omission_FPR,
                            no_omission_MCR, no_omission_NPP, no_omission_ODP,
                            no_omission_OR, no_omission_PPP, no_omission_prevalence,
                            no_omission_TNR, no_omission_TPR)
  
  # Create a dataframe of evaluation statistics using the threshold at which 
  # the modeled prevalence is closest to the observed prevalence
  prevalence_prevalence = t(as.data.frame(as.double(slot(evaluation_prevalence,"prevalence"))))
  prevalence_OR = t(as.data.frame(as.double(slot(evaluation_prevalence,"OR"))))
  prevalence_MCR = t(as.data.frame(as.double(slot(evaluation_prevalence,"MCR"))))
  prevalence_NPP = t(as.data.frame(as.double(slot(evaluation_prevalence,"NPP"))))
  prevalence_PPP = t(as.data.frame(as.double(slot(evaluation_prevalence,"PPP"))))
  prevalence_FNR = t(as.data.frame(as.double(slot(evaluation_prevalence,"FNR"))))
  prevalence_FPR = t(as.data.frame(as.double(slot(evaluation_prevalence,"FPR"))))
  prevalence_TNR = t(as.data.frame(as.double(slot(evaluation_prevalence,"TNR"))))
  prevalence_TPR = t(as.data.frame(as.double(slot(evaluation_prevalence,"TPR"))))
  prevalence_CCR = t(as.data.frame(as.double(slot(evaluation_prevalence,"CCR"))))
  prevalence_ODP = t(as.data.frame(as.double(slot(evaluation_prevalence,"ODP"))))
  
  prevalence_bind <- rbind(prevalence_CCR, prevalence_FNR, prevalence_FPR,
                           prevalence_MCR, prevalence_NPP, prevalence_ODP,
                           prevalence_OR, prevalence_PPP, prevalence_prevalence,
                           prevalence_TNR, prevalence_TPR)
  
  # Create a dataframe of evaluation statistics using the threshold where the
  # sensitivity and specificity are equal
  equal_sens_spec_prevalence = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"prevalence"))))
  equal_sens_spec_OR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"OR"))))
  equal_sens_spec_MCR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"MCR"))))
  equal_sens_spec_NPP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"NPP"))))
  equal_sens_spec_PPP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"PPP"))))
  equal_sens_spec_FNR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"FNR"))))
  equal_sens_spec_FPR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"FPR"))))
  equal_sens_spec_TNR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"TNR"))))
  equal_sens_spec_TPR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"TPR"))))
  equal_sens_spec_CCR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"CCR"))))
  equal_sens_spec_ODP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"ODP"))))
  
  equal_sens_spec_bind <- rbind(equal_sens_spec_CCR, equal_sens_spec_FNR,
                                equal_sens_spec_FPR, equal_sens_spec_MCR, equal_sens_spec_NPP,
                                equal_sens_spec_ODP, equal_sens_spec_OR, equal_sens_spec_PPP,
                                equal_sens_spec_prevalence, equal_sens_spec_TNR, equal_sens_spec_TPR)
  
  # Generate evaluation statistics table
  all_thresholds_matrix <- cbind(kappa_bind, spec_sens_bind, no_omission_bind,
                                 prevalence_bind, equal_sens_spec_bind)
  
  # Add in the original threshold values
  colnames(all_thresholds_matrix) <- colnames(threshold_all)
  threshold_matrix <- rbind(threshold_all, all_thresholds_matrix) 
  
  # Change column and row names
  colnames(threshold_matrix) <- c("Max Kappa", "Max sum sensitivity and specificity",
                                  "No Omission", "Prevalence", "Equal sensitivity specificity") 
  rownames(threshold_matrix) <- c("Threshold value", "Prevalence","Odds ratio",
                                  "Misclassification rate","Negative predictive power",
                                  "Positive predictive power","False negative rate",
                                  "False positive rate","True negative rate",
                                  "True positive rate","Correct classification rate",
                                  "Overall diagnostic power")
  
  # Write evaluation statistics as csv
  write.csv(threshold_matrix, file=file.path(EC.env$outputdir, paste0('evaluation_',
                                                                      species, '_maxent.csv')))
  
  # return equal sensitivity specificity threshold
  return(threshold_equal_sens_spec)
}
