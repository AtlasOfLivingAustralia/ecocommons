#' read_json
#' 
#' Local function for importing json files. Doesn't do much
#' 
#' @param file a file path
#'
#' @export EC_read_json

EC_read_json <- function(file){
  
  result <- rjson::fromJSON(file = file)
  class(result) <- "model_parameters"
  return(result)
  
}


print.model_parameters <- function(x){
  str(x, max.level = 2)
}

