#' Function to save the  graphs generated between multiple devices?
#'
#' @param fileroot 
#' @param ext 
#'
#' @export EC_DevSave
#' @importFrom dev2 dev.copy2eps
#' @importFrom dev2 dev.copy2pdf
#'

EC_DevSave <- function(fileroot, ext=".pdf") {
  if (ext==".eps") {dev2::dev.copy2eps(file=paste(fileroot,ext,sep="."))
  
    } else {
      dev2::dev.copy2pdf(file=paste(fileroot,"pdf",sep="."))}
}
