#' Retrieve coordinate projection from both species user and climate data
#' 
#' This script compares the projection from the species user data with the
#' climate data projections, transforming as saving to the same projection if
#' necessary, and creating an outfile name
#' 
#' @param data 
#' @param climate.data 
#'
#' @importFrom raster compareCRS
#' @importFrom sp SpatialPoints
#' @importFrom sp spTransform
#' 
#' @export EC_data_projection

EC_data_projection <- function(data,
                               climate.data) {
  
  if (!is.null(data) & !raster::compareCRS(data, climate.data, verbatim=TRUE)) {
    species <- sp::SpatialPoints(data)
    if (is.na(crs(species))) {
      crs(species) <- '+init=epsg:4326'
    }
    newdata <- as.data.frame(sp::spTransform(species, crs(climate.data))) # convert to data frame
    names(newdata) <- names(data)
    return(newdata)
  }
  return(data)
}
