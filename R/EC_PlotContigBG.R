#' Load and export plot function for continuous predictor variables for the 
#' background data
#'
#' @param dataToExplore 
#' @param listVars 
#' @param doCorrelation 
#' @param outerTitle 
#' @param fnamePrefix 
#' @param outputdir 
#'
#' @export EC_PlotContinBG
#' @importFrom 
#' 

EC_PlotContinBG <- function(dataToExplore, listVars, doCorrelation=TRUE, 
                                outerTitle="", fnamePrefix="background_predictors_", outputdir="output")  {
  nVars = length(listVars)
  
  if (nVars == 0) {
    return()
  }
  
  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)
  plotWidth = 800; plotHeight = 800
  
  # Boxplots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'boxplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    boxplot(dataToExplore[,c(nameVar)], main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Histograms
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'histogram.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    hist(dataToExplore[,c(nameVar)], breaks=10, xlab="", main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Density plots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'density.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    plot(density(dataToExplore[,c(nameVar)], na.rm=TRUE), main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
}
