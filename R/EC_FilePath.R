# function to generate a filename for the specified file type and extension.
EC_FilePath <- function(file_type, projection_name, species, 
                        outputdir=EC.env$outputdir, filename_ext=NULL, file_ext='tif') {
  if (is.null(filename_ext)) {
    basename = paste(file_type, projection_name, species, sep="_")
  }  else {
    basename = paste(file_type, projection_name, species, filename_ext, sep="_")
  }
  return(file.path(outputdir, paste(basename, file_ext, sep=".")))
}