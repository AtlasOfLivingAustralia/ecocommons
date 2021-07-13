#' Load and export plot function for categorical predictor variables for
#' the occurrence data
#' 
#' @param dataToExplore 
#' @param listVars 
#' @param outerTitle 
#' @param fnamePrefix 
#' @param outputdir 
#'
#' @export EC_PlotCategOCC
#' 
#' 

EC_PlotCategOCC <- function(dataToExplore, listVars, outerTitle="", 
                            fnamePrefix="occurrence_predictors_", outputdir="output") {
  
  nVars = length(listVars)
  
  if (nVars == 0) {
    return()
  }
  
  plotWidth = 800; plotHeight = 800
  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)
  
  # Barplots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'barplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    barplot(table(dataToExplore[,c(nameVar)]), 
            main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
}