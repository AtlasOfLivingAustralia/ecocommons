#' Function to crop the raster to the extent of the constraint region and 
#' mask it
#' 
#' @param raster 
#' @param parsed_geojson 
#'

#' @export
#' @importFrom raster crop
#' @importFrom raster mask
#'

EC_Mask <- function(raster, parsed_geojson) {
  
  cropped_raster <- raster::crop(raster, extent(parsed_geojson))  # crop to constraint region before masking
  
  envraster_filename <- paste(EC.env$workdir, basename(tempfile(fileext = ".grd")), sep="/")
  masked_raster <- raster::mask(cropped_raster, parsed_geojson, filename = envraster_filename)
  
  if (is.factor(masked_raster)) {
    masked_raster = as.factor(masked_raster)
  }
  
  EC_RevRasterObject(stack(cropped_raster))  # remove cropped raster and associated raster files (i.e. grd and gri)
  
  return(masked_raster)
}