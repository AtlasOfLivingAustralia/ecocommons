# load raster and assign projection (crs) if missing

EC_RasterLoad <- function(filename) {
  r = raster(filename)
  if (is.na(crs(r))) {
    crs(r) = CRS("+init=epsg:4326")
  }
  return(r)
}