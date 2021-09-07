#' Function to generate a file name for the specified file type and extension
#'
#' @param file_type 
#' @param projection_name 
#' @param species 
#' @param outputdir 
#' @param filename_ext 
#' @param file_ext 

EC_file_path <- function (file_type, projection_name, species, 
                          outputdir = EC.env$outputdir, filename_ext=NULL, 
                          file_ext='tif') {
  
  if (is.null(filename_ext)) {
    basename = paste(file_type, projection_name, species, sep="_")
  }  else {
    basename = paste(file_type, projection_name, species, filename_ext, sep="_")
  }
  
  return (file.path(outputdir, paste(basename, file_ext, sep=".")))
}