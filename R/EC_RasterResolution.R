# Function to determine and resample the resolution of the user raster
#
# resamplingflag: a flag to determine which resampling approach to take
# selected_layers: a list of indexes to the raster layers to be considered when determine the resolution to be used.
EC_RasterResolution <- function(rasters, resamplingflag, selected_layers) {
  get.resolution <- function(i, rasters) {
    # Return the resolution of the raster given by the index
    return(res(rasters[[i]]))
  }
  
  if (is.null(selected_layers)) {
    resolutions = lapply(rasters, res) # Get resolutions of the raster layers
  }
  else {
    resolutions = lapply(selected_layers, get.resolution, rasters) # get the resolutions of the of the selected raster layers only 
  }
  
  if (resamplingflag == "highest") {
    common.res = Reduce(pmin, resolutions)
  } else if (resamplingflag == "lowest") {
    common.res = Reduce(pmax, resolutions)
  }
  
  resolutions = lapply(rasters, res) # Get resolutions of all input raster layers
  is.same.res = all(sapply(resolutions, function(x) all(common.res == x)))
  return (list(common.res=common.res, is.same.res=is.same.res))
}