#' Function to save RData in outputdir
#'
#' @param robj 
#' @param name 
#' @param outputdir 
#'
#' @export EC_save
#' 
#' 

EC_save <- function(robj, name, outputdir=EC.env$outputdir) {
  filename = file.path(outputdir, name)
  save(robj, file=filename)
}