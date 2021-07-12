#' Apply same extension and projection for all rasters
#' rasters: a vector of rasters, all rasters should have same resolution
#' common.crs: crs to use to calculate intersection
#' @param rasters 
#' @param common.crs 
#'
#' @export
#' @importFrom raster intersect
#'

EC_RasterExtent <- function(rasters, common.crs) {
  extent.list = lapply(rasters, function(r) { extent(projectExtent(r, common.crs)) }) # intersect all extents
  
  common.extent = Reduce(raster::intersect, extent.list) # if all extents are the same (used to print warning)
  
  equal.extents = all(sapply(extent.list, function (x) common.extent == x))
  
  return (list(equal.extents = equal.extents, common.extent = common.extent))
}


