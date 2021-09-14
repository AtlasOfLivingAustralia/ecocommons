#' Load raster and assign projection (crs) if missing
#'
#' @param filename 
#'
#' @importFrom raster raster
#' @importFrom raster crs
#' 
#' @export EC_read_raster

EC_read_raster <- function(filename) {
  r = raster::raster(filename)
  if (is.na(crs(r))) {
    raster::crs(r) <- crs("+init=epsg:4326")
  }
  return(r)
}
