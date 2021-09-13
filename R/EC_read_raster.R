#' Load raster and assign projection (crs) if missing
#'
#' @param filename 
#'
#' @importFrom raster raster
#' @importFrom sp CRS
#' 
#' @export EC_read_raster

EC_read_raster <- function(filename) {
  r = raster::raster(filename)
  if (is.na(crs(r))) {
    crs(r) <- sp::CRS("+init=epsg:4326")
  }
  return(r)
}
