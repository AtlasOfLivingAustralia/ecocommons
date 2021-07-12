# Project species distribution data and climate data in the same crs (EPSG:4326)

EC_SpDataProjection <- function(data, climate.data) {
  sp <- SpatialPoints(data)
  if (is.na(crs(sp))) {
    crs(sp) <- '+init=epsg:4326'
  }
  
  if (!compareCRS(sp, climate.data, verbatim=TRUE)) {
    sp <- spTransform(sp, crs(climate.data))
  }
  return(sp)
}