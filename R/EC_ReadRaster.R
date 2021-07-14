#' Load raster and assign projection (crs) if missing
#'
#' @param filename 
#'
#' @export EC_ReadRaster
#' @importFrom raster raster
#' @importFrom sp CRS
#' 
#' 

EC_ReadRaster <- function(filename) {
  r = raster(filename)
  if (is.na(crs(r))) {
    crs(r) = CRS("+init=epsg:4326")
  }
  return(r)
}
