#' Save pdf
#'
#' @param ... 
#' @param filename 
#' @param aspdf 
#' @param outputdir 
#'
#' @importFrom gridExtra grid.arrange
#'
#' @export EC_save_pdf
#'

EC_save_pdf <- function(..., 
                        filename,
                        aspdf,
                        outputdir = EC.env$outputdir) {
  # Save output in PDF format 
  
  if (aspdf) {
    png(file= file.path(outputdir, paste(filename, 'png', sep = ".")))
    gridExtra::grid.arrange(...)
    dev.off()
  } else {
    gridExtra::grid.arrange(...)
  }
}
