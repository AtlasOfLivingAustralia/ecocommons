# function to save projection output raster

EC_SaveModelProjection <- function(model.obj, projection.name, species, 
                                     species_algo_str, outputdir=EC.env$outputdir, filename_ext=NULL) {
  
  filename = EC_FilePath("proj", projection.name, species_algo_str, outputdir, filename_ext, "tif")
  writeRaster(model.obj, filename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  pngfilename = EC_FilePath("proj", projection.name, species_algo_str, outputdir, filename_ext, "png")
  png(pngfilename)
  title = paste(species, projection.name, "projections", sep=" ")
  plot(model.obj, xlab="latitude", ylab="longtitude", main=title)
  dev.off()
  
  return (filename)
}