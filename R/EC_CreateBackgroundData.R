#' Generate background data on the region of interest
#'
#' @param rastermask 
#' @param bgsize 
#'
#' @export EC_CreateBackgroundData
#' @ImportFrom dismo randomPoints
#'

EC_CreateBackgroundData <- function(rastermask, bgsize) {
  background <- dismo::randomPoints(rastermask, n=bgsize, extf=1.0, excludep=FALSE)
  colnames(background) <- c('lon', 'lat')
  return(as.data.frame(background))
}