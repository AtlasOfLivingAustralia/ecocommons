# function to save CSV Data in outputdir
EC_WriteCSV <- function(robj, name, outputdir=EC.env$outputdir, rownames=TRUE) {
  filename = file.path(outputdir, name)
  write.csv(robj, file=filename, row.names=rownames)
}