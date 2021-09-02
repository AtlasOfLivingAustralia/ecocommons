#' SAve output in PDF format
#'
#' @param ... 
#' @param filename 
#' @param aspdf 
#' @param outputdir 
#'
#' @export EC_save_pdf
#' @importFrom gridExtra grid.arrange

EC_save_pdf <- function(..., filename, aspdf, outputdir=EC.env$outputdir)
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
