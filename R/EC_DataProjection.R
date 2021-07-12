# Retrieve coordinate projection from both species user and climate data

# This script compares the projection from the species user data with the
# climate data projections, transforming as saving to the same projection if 
# necessary, and creating an outfile name

EC_DataProjection <- function(data, climate.data) {
  if (!is.null(data) & !compareCRS(data, climate.data, verbatim=TRUE)) {
    sp <- SpatialPoints(data)
    if (is.na(crs(sp))) {
      crs(sp) <- '+init=epsg:4326'
    }
    newdata <- as.data.frame(spTransform(sp, crs(climate.data))) # convert to data frame
    names(newdata) <- names(data)
    return(newdata)
  }
  return(data)
}

EC_OutfileName <- function(filename, id_str, ext) {
  return(sprintf("%s_%s.%s", filename, id_str, ext))
}
