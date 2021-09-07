#' Function to save projection as png image
#'
#' @param inputfile 
#' @param projection_name 
#' @param species 
#' @param species_algo_str 
#' @param outputdir 
#' @param filename_ext 
#' 
#' @importFrom raster raster

EC_save_projection_image <- function(inputfile, projection_name, species, 
                                     species_algo_str,
                                     outputdir = EC.env$outputdir,
                                     filename_ext = NULL) {
  
  filename <- EC_file_path ("proj", projection_name, species_algo_str, 
                            outputdir, filename_ext, "png")
  png(filename)
  
  title = paste(species, projection_name, "projections", sep=" ")
  
  plot(raster::raster(inputfile), main = title,
       xlab = 'longitude', ylab = 'latitude')
  # TODO: to use levelplot to produce histogram instead of plot.
  #EC_plot_projection(inputfile, title)
  dev.off()
}