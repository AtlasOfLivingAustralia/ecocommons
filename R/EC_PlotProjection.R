#' Function to plot projection tiff file (with histogram)
#'
#' @param inputfile 
#' @param main 
#'
#' @importFrom rasterVis levelplot
#' 
#' 

EC_PlotProjection <- function(inputfile, main) {
  my.at <- seq(0,1.0, by=0.1)  # breaks of the color key
  my.labs.at <- seq(0,1.0, by=0.25)
  my.lab <- seq(0,1.0, by=0.25)
  my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)
  
  # read in tiff input file as rasterstack and plot it
  require('rasterVis')
  levelplot(stack(raster(inputfile)),
            at=my.at,
            margin=T,
            col.regions=my.col,
            main=main,
            colorkey=list(labels=list(
              labels=my.lab,
              at=my.labs.at)))
}