#' Function to save the  graphs generated between multiple devices?
#'
#' @param fileroot 
#' @param ext 
#'
#' @export EC_save_dev

EC_save_dev <- function(fileroot, ext=".pdf") {
  
  if (ext==".eps") {dev.copy2eps(file=paste(fileroot,ext,sep="."))
    } else {
      dev.copy2pdf(file=paste(fileroot,"pdf",sep="."))}
}
