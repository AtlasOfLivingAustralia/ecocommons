#' Function to save CSV Data in outputdir
#'
#' @param robj 
#' @param name 
#' @param outputdir 
#' @param rownames 
#'
#' @export EC_write_csv

EC_write_csv <- function(robj, name, outputdir=EC.env$outputdir, rownames=TRUE) {
  
  filename = file.path(outputdir, name)
  
  write.csv(robj, file=filename, row.names=rownames)
}