#' Function to get model object
#'
#' @param model.file 
#'
#' @export
#' @importFrom 
#' 

EC_GetModelObject <- function(model.file=EC.env$inputmodel) {
  return (get(load(file=model.file)))
}
