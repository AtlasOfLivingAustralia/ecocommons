#' Function to compute and save occurrence probability change metrics as geotif file
#'
#' @param prob_rasters 
#' @param outfilename 
#'
#' @export EC_OccurenceProbMetric
#' @importFrom raster overlay
#' @importFrom raster writeRaster
#' 

EC_OccurenceProbMetric <- function(prob_rasters, outfilename) {
  
  changeproj <- raster::overlay(prob_rasters[[1]], prob_rasters[[2]], 
                        fun=function(r1, r2) { return(r1-r2) })
  
  raster::writeRaster(changeproj, outfilename, format="GTiff", 
              options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
}
