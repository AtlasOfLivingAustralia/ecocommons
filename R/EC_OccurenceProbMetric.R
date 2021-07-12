# function to compute and save occurrence probability change metrics as geotif file
EC_OccurenceProbMetric <- function(prob_rasters, outfilename) {
  changeproj <- overlay(prob_rasters[[1]], prob_rasters[[2]], 
                        fun=function(r1, r2) { return(r1-r2) })
  writeRaster(changeproj, outfilename, format="GTiff", 
              options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
}
