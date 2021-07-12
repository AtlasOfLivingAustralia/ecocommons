# function to get model object
EC_GetModelObject <- function(model.file=EC.env$inputmodel) {
  return (get(load(file=model.file)))
}
