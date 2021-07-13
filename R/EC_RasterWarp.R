#' This function just need to be used if raster layers are not ligned up
#' correctly.
#'
#' @param raster.filenames 
#' @param raster.types 
#' @param reference 
#' @param overwrite 
#'
#' @export
#' @importFrom gdalUtils gdalwarp
#' @importFrom raster raster
#' @importFrom rgdal GDALinfo
#' 

EC_RasterWarp <- function(raster.filenames, raster.types, reference, overwrite=TRUE) {

   rasters <- mapply(
    function(filename, filetype) {
      gdinfo <- rgdal::GDALinfo(filename)
      mdata <- attr(gdinfo, 'df')
      dtype <- as.character(mdata[['GDType']])
      hasNoDataValues <- mdata[['hasNoDataValue']]
      
      r <- EC_RasterLoad(filename)  # warp, crop and rescale raster file if necessary
      dir <- dirname(filename)
      temp_raster <- file.path(dir, paste0(basename(tempfile()), '.tif'))
      te <- extent(reference)
      
      # This is to fix issue with NA value being treated as value 0 if nodatavalue is not set.
      if (hasNoDataValues) {
        gdalUtils::gdalwarp(filename, temp_raster,
                 s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                 te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                 ts=c(ncol(reference), nrow(reference)),
                 # tr=c(...), ... either this or ts
                 r="near",
                 of="GTiff",
                 dstnodata=mdata[['NoDataValue']],
                 co=c("TILED=YES", "COMPRESS=LZW")
        )
      }
      else {
        gdalUtils::gdalwarp(filename, temp_raster,
                 s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                 te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                 ts=c(ncol(reference), nrow(reference)),
                 # tr=c(...), ... either this or ts
                 r="near",
                 of="GTiff",
                 co=c("TILED=YES", "COMPRESS=LZW")
        )
      }
      
      rasterfilename <- temp_raster
      if (overwrite) {
        file.rename(temp_raster, filename)
        rasterfilename = filename 
      }

      r <- raster::raster(rasterfilename)
      if (filetype == "categorical") { # convert to factor if categorical
        r = as.factor(r)
      }
      return(r)
    },
    raster.filenames, raster.types)
  return(rasters)
}