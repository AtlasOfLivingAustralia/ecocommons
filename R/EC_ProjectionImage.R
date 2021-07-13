#' Function to save projection as png image
#'
#' @param inputfile 
#' @param projection.name 
#' @param species 
#' @param species_algo_str 
#' @param outputdir 
#' @param filename_ext 
#'
#' 

EC_ProjectionImage <- function(inputfile, projection.name, species, 
                                     species_algo_str, outputdir=EC.env$outputdir, filename_ext=NULL) {
  
  filename <- EC_FilePath("proj", projection.name, species_algo_str, 
                         outputdir, filename_ext, "png")
  png(filename)
  
  title = paste(species, projection.name, "projections", sep=" ")
  
  plot(raster(inputfile), main=title, xlab='longitude', ylab='latitude')
  # TODO: to use levelplot to produce histogram instead of plot.
  #EC_PlotProjection(inputfile, title)
  dev.off()
}