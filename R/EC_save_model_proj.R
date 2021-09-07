#' Function to save projection output raster
#'
#' @param model_obj 
#' @param projection_name 
#' @param species 
#' @param species_algo_str 
#' @param outputdir 
#' @param filename_ext 
#'
#' @importFrom raster writeRaster
#'
#' @export EC_save_model_proj

EC_save_model_proj <- function(model_obj, projection_name, species,
                               species_algo_str, outputdir = EC.env$outputdir, 
                               filename_ext = NULL) {
  
  filename = EC_file_path("proj", projection_name, species_algo_str, outputdir, 
                          filename_ext, "tif")
  
  raster::writeRaster(model_obj, filename, format="GTiff",
              options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  
  pngfilename = EC_file_path("proj", projection_name, species_algo_str, 
                             outputdir, filename_ext, "png")
  
  png(pngfilename)
  
  title = paste(species, projection_name, "projections", sep=" ")
  
  plot(model_obj, xlab="latitude", ylab="longtitude", main=title)
  
  dev.off()
  
  return (filename)
}