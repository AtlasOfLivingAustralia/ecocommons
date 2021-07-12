# Stack raster input files from user environment
# selected_layers is a list of layers to be considered when determine the
# resolution of the raster. If none, consider all layers.

EC_EnviroStack <- function(filenames, types, layernames, resamplingflag,
                           selected_layers=NULL) {
  rasters <- EC_RasterResampled (filenames, types, resamplingflag,
                                 selected_layers) # adjust rasters to same projection, resolution and extent

  rasterstack <- stack(rasters)

  names(rasterstack) <- unlist(layernames)  # assign predefined variable names
  return(rasterstack)
}
