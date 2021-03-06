#' Function to save projection 
#'
#' @param proj.model 
#' @param species_algo_str 
#' @param filename_ext 
#'
#' @export EC_save_projection

EC_save_projection <- function(proj.model,
                               species_algo_str,
                               filename_ext=NULL) {
  
  if (!is.null(filename_ext)) {
    basename = paste ("proj", 'current', species_algo_str, filename_ext, sep="_")
    
    } else {
      
      basename = paste ("proj", 'current', species_algo_str, sep="_")
  }
  
  png(file=file.path(EC.env$outputdir, paste(basename, 'png', sep=".")))
  
  plot(proj.model, on_0_1000=FALSE)
  
  dev.off()
}
