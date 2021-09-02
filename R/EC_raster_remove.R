#' Remove raster object and its associated raster files (i.e. grd and gri) if any
#'
#' #@param rasterObject 
#'
#' #@importFrom spatial.tools raster_to_filenames

EC_raster_remove <- function(rasterObject) {
  raster_filenames <- raster::filename(rasterObject)
  for (fname in raster_filenames) {
    if (extension(fname)  == '.grd') {
      file.remove(fname, extension(fname, '.gri'))
    }
  }
  rm(rasterObject)
}
