#' This function aims to create the same projection for all the layers.
#' If a common projection is not established, it will be projected in EPSG:4326,
#' common in EcoCommons
#'
#' @param rasters 
#' @param resamplingflag 
#' @param selected_layers 
#'
#' @export
#' @importFrom raster crs
#' @importFrom sp CRS
#' 
#' 

EC_RasterRef <- function(rasters, resamplingflag, selected_layers) {
  empty.rasters <- lapply(rasters, function(x) { projectExtent(x, crs(x)) })  # create list of empty rasters
  common.crs <- crs(empty.rasters[[1]]) # choose a common.crs if all crs in rasters are the same use that one

  if (! do.call(compareRaster, c(empty.rasters, extent=FALSE, rowcol=FALSE, 
                                 prj=TRUE, res=FALSE, orig=FALSE, rotation=FALSE, stopiffalse=FALSE))) {
    common.crs = CRS("+init=epsg:4326") # common projection adopted in EcoCommons
    EC_LogWarning(sprintf("Auto projecting to common CRS: %s", common.crs))
    empty.rasters = lapply(empty.rasters, function(x) { projectExtent(x, common.crs) })
  }

  ce <- EC_RasterExtent(empty.rasters, common.crs)
  if (! ce$equal.extents) {
    EC_LogWarning(sprintf("Auto cropping to common extent %s", EC_RasterExtentToSTR(ce$common.extent)))
  }

  cr <- EC_RasterResolution(empty.rasters, resamplingflag, selected_layers)
  if (! cr$is.same.res) {
    EC_LogWarning(sprintf("Auto resampling to %s resolution [%f %f]", resamplingflag,
                          cr$common.res[[1]], cr$common.res[[2]]))
  }

  empty.rasters <- lapply(
    empty.rasters,
    function(x) {
      extent(x) = ce$common.extent
      res(x) = cr$common.res
      return(x)
    })

  return(empty.rasters[[1]])  # from now on all empty.rasters should be in the same projection
}
