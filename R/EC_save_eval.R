#' Save a evalution model for all statistics
#'
#' @param out_evaluation 
#' @param out_stats 
#' @param out_lossfunction 
#' @param species_algo_str 
#'
#' @export EC_save_eval

EC_save_eval <- function(out_evaluation, out_stats, out_lossfunction, 
                         species_algo_str){
  
  EC_write_csv(data.frame(out_evaluation),
              name = paste0(paste("Evaluation-data", species_algo_str, 
                                  sep="_"), ".csv"))
  
  EC_write_csv(data.frame(out_stats),
              name = paste0(paste("Evaluation-statistics", species_algo_str, 
                                  sep="_"), ".csv"))
  
  EC_write_csv(data.frame(out_lossfunction),
              name = paste0(paste("Loss-function-intervals-table", 
                                  species_algo_str, sep="_"), ".csv"))
}
