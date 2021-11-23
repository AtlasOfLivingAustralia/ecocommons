#' Function to save projection 
#'
#' @param model_proj 
#' @param species_algo_str 
#' @param filename_ext 
#' 
#' @importFrom raster plot
#'
#' @export EC_save_projection

EC_save_projection <- function(model_proj,
                               species_algo_str,
                               filename_ext = NULL) {
  
  if (!is.null(filename_ext)) {
    basename = paste ("proj", 'current', species_algo_str, filename_ext, sep = "_")
    
    } else {
      
      basename = paste ("proj", 'current', species_algo_str, sep = "_")
  }
  
  png(file = file.path(EC.env$outputdir, paste(basename, 'png', sep = ".")))
  
  raster::plot (model_proj)
  
  dev.off()
}
