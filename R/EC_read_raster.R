#' Load raster and assign projection (crs) if missing
#'
#' @param filename 
#'
#' @export EC_read_raster
#' @importFrom raster raster
#' @importFrom sp CRS

EC_read_raster <- function(filename) {
  r = raster::raster(filename)
  if (is.na(crs(r))) {
    crs(r) <- CRS("+init=epsg:4326")
  }
  return(r)
}
