#' If chosen to use a bias file, apply it in the occurences to the background
#' points
#' 
#' @param bg 
#' @param biasfile 
#' @param bgsize 

#' @export EC_CreateBiasedBackgroundData
#' @importFrom raster extract
#' @importFrom raster raster
#' @importFrom sp SpatialPoints

EC_CreateBiasedBackgroundData <- function(bg, biasfile, bgsize) {
  
  bias <- raster::raster(biasfile)

  bg_sp_class<-sp::SpatialPoints(coords=cbind(bg$lon,bg$lat),
                                 proj4string=CRS(as.character(NA)))   # convert background points to spatial points formal class
  
  weight <- raster::extract(bias, bg_sp_class)  
  bg_sp_class <- bg_sp_class[!is.na(weight), ]
  weight <- weight[!is.na(weight)]  # random sample of the candidate background points
  bg_biased <- bg_sp_class[sample(length(bg_sp_class), size=bgsize, replace = TRUE, prob=weight), ]
  colnames(bg_biased@coords) <- c("lon", "lat")
  bg_biased <- bg_biased@coords
  return(bg_biased)
}
