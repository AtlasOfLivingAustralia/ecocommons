# function to save RData in outputdir

EC_Save <- function(robj, name, outputdir=EC.env$outputdir) {
  filename = file.path(outputdir, name)
  save(robj, file=filename)
}