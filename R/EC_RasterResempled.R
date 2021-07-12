# Load rasters, determine common projection and convert categorical
# rasters to factors
# components:
# raster.filenames is a vector of filenames that will be loaded as rasters
# resamplingflag is a flag to determine which resampling approach to take

EC_RasterResampled  <- function(raster.filenames, raster.types, resamplingflag, selected_layers=NULL, overwrite=TRUE) {
  rasters <- lapply(raster.filenames, EC_RasterLoad)
  
  reference <- EC_RasterRef(rasters, resamplingflag, selected_layers) # determine common raster shape
  
  rasters <- EC_RasterWarp(raster.filenames, raster.types, reference, overwrite) # adjust rasters spatially and convert categorical rasters to factors
  
  return(rasters)
}