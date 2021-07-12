EC_SaveModelEval <- function(out.evaluation, out.stats, out.lossfunction, species_algo_str){
  EC_WriteCSV(data.frame(out.evaluation),
              name = paste0(paste("Evaluation-data", species_algo_str, sep="_"), ".csv"))
  EC_WriteCSV(data.frame(out.stats),
              name = paste0(paste("Evaluation-statistics", species_algo_str, sep="_"), ".csv"))
  EC_WriteCSV(data.frame(out.lossfunction),
              name = paste0(paste("Loss-function-intervals-table", species_algo_str, sep="_"), ".csv"))
}
