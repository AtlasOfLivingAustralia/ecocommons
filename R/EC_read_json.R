#' Local function for importing json files
#' 
#' @param file a file path
#' 
#' @importFrom rjson fromJSON
#' 
#' @export EC_read_json

EC_read_json <- function (file){
  
  result <- rjson::fromJSON (file = file)
  class(result) <- "model_parameters"
  return(result)
  
}


#_____________________________________________________________________________


print.model_parameters <- function (x){
  str(x, max.level = 2)
}

