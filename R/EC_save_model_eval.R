#' Save evaluation models for all statistics in csv format
#'
#' @param out.evaluation 
#' @param out.stats 
#' @param out.lossfunction 
#' @param species_algo_str 
#'
#' @export EC_save_model_eval

EC_save_model_eval <- function(out.evaluation,
                               out.stats,
                               out.lossfunction,
                               species_algo_str){
  
  EC_write_csv(data.frame(out.evaluation),
              name = paste0(paste("Evaluation-data",
                                  species_algo_str, sep="_"), ".csv"))
  
  EC_write_csv(data.frame(out.stats),
              name = paste0(paste("Evaluation-statistics",
                                  species_algo_str, sep="_"), ".csv"))
  
  EC_write_csv(data.frame(out.lossfunction),
              name = paste0(paste("Loss-function-intervals-table",
                                  species_algo_str, sep="_"), ".csv"))
}
