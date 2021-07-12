# rasters: a vector of rasters, all rasters should have same resolution
# common.crs: crs to use to calculate intersection
EC_RasterExtent <- function(rasters, common.crs) {
  extent.list = lapply(rasters, function(r) { extent(projectExtent(r, common.crs)) }) # intersect all extents
  common.extent = Reduce(raster::intersect, extent.list) # compare all against common extents to find out if all extents are the same (used to print warning)
  equal.extents = all(sapply(extent.list, function (x) common.extent == x))
  return (list(equal.extents=equal.extents, common.extent=common.extent))
}

EC_RasterExtentToSTR <- function(ext) {
  return(sprintf("xmin=%f xmax=%f ymin=%f ymax=%f", ext@xmin, ext@xmax, ext@ymin, ext@ymax));
}
