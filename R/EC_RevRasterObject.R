#' Remove raster object and its associated raster files (i.e. grd and gri) if any
#'
#' #@param rasterObject 
#'
#' #@importFrom spatial.tools raster_to_filenames
#' # raster filename ### OR ### terra sources

EC_RevRasterObject <- function(rasterObject) {
  raster_filenames <- raster_to_filenames(rasterObject, unique = TRUE)
  for (fname in raster_filenames) {
    if (extension(fname)  == '.grd') {
      file.remove(fname, extension(fname, '.gri'))
    }
  }
  rm(rasterObject)
}