#' SAve output in PDF format
#'
#' @param ... 
#' @param filename 
#' @param aspdf 
#' @param outputdir 
#'
#' @export EC_SavePDF
#' @importFrom gridExtra grid.arrange
#'

EC_SavePDF <- function(..., filename, aspdf, outputdir=EC.env$outputdir)
{
  library("gridExtra")
  if (aspdf)
  {
    png(file=file.path(outputdir, paste(filename, 'png', sep=".")))
    grid.arrange(...)
    dev.off()
  }
  else {
    grid.arrange(...)
  }
}
