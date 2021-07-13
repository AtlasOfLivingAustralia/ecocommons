#' This is a fix for SDM tool distance function due to floating point operation
#'
#' @param lat1 
#' @param lon1 
#' @param lat2 
#' @param lon2 
#' @param bearing 
#'
#' @export EC_Distance
#' @importFrom raster atan2
#' @importFrom SDMTools
#' 

EC_Distance <- function(lat1, lon1=NULL, lat2=NULL, lon2=NULL, bearing=FALSE) {
  if (is.data.frame(lat1) | is.matrix(lat1)) {  # if yes, break it out to individual vectors
    lat1 = as.matrix(lat1); if (ncol(lat1)!=4) 
      stop('incorrect lat/lon inputs... must be matrix with 4 columns or 4 vectors')
    lon2=lat1 [,4]; lat2=lat1 [,3]; lon1=lat1 [,2]; lat1=lat1 [,1]  # break out individual columns
    
  } else if (!is.null(lat2) & !is.null(lon1) & !is.null(lon2)) {
    if (!all(c(length(lat2),length(lon1),length(lon2))==length(lat1)))
      stop('inputs must all be of same length')
    
  } else { stop('inappropriate inputs... see helpfile') }
  
  if (any(c(lon1,lon2) < - 180.0001) | any(c(lon1,lon2) > 180.0001))
    stop('lon must be decimal degrees between -180 & 180')
  if (any(c(lat1,lat2) < - 90.0001) | any(c(lat1,lat2) > 90.0001))
    stop('lat must be decimal degrees between -90 & 90')
  
  # Cycle through and output the new data
  output <- data.frame(lat1=lat1, lon1=lon1, lat2=lat2, lon2=lon2)
  output$distance <- round(.Call('Dist',output$lat1, output$lon1, output$lat2,
                                 output$lon2, PACKAGE='SDMTools'),2)  # round to the nearest mm
  if (bearing) {  # if requested, calculate bearing
    lat1<- lat1 * pi/180; lat2= lat2 * pi/180; lon1=lon1 * pi/180;
    lon2= lon2 * pi/180  # convert to radians
    brng <-  atan2(sin(lon2 - lon1) * cos(lat2),cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2))  # estimate bearing
    output$bearing = ((brng * 180/pi) + 360) %%360  # convert to bearing in degrees
  }

  return(output)  # return the output
}