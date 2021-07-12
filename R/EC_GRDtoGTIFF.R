# convert all .gri/.grd found in folder to gtiff

EC_GRDtoGTIFF <- function(folder, algorithm, filename_ext=NULL, noDataValue=NULL) {
  grdfiles <- list.files(path=folder,
                         pattern="^.*\\.grd")
  for (grdfile in grdfiles) {
    
    ext = file_ext(grdfile)
    if (!is.null(ext)) {
      pattern = paste0('\\.', ext, '$')
      grdname <- sub(pattern, '', grdfile)
    }
    
    grd <- raster(file.path(folder, grdfile))
    
    if (is.na(proj4string(grd))) {
      crs = CRS("+init=epsg:4326")  # if projection is missing, initialise it to EPSG:4326
      proj4string(grd) <- crs
    }
    
    basename = paste(grdname, algorithm, sep="_")
    if (!is.null(filename_ext)) {
      basename = paste(grdname, algorithm, filename_ext, sep="_")
    }
    filename = file.path(folder, paste(basename, 'tif', sep="."))
    
    dtype = dataType(grd)
    if (is.null(noDataValue)) {
      writeRaster(grd, filename, datatype=dataType(grd),
                  format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
    }
    else {
      writeRaster(grd, filename, datatype=dataType(grd), NAflag=noDataValue,
                  format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
    }
    file.remove(file.path(folder, paste(grdname, c("grd","gri"), sep=".")))  # remove grd files
  }
}
