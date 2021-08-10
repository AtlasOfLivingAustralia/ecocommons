#' Stack raster input files from user environment
#' selected_layers is a list of layers to be considered when determine the
#' resolution of the raster. If none, consider all layers.
#'
#' @param filenames 
#' @param types 
#' @param layernames 
#' @param resamplingflag 
#' @param selected_layers 
#'
#' @export EC_EnviroStack
#' @importFrom raster compareRaster
#' @importFrom raster stack
#' @importFrom raster crs
#' @importFrom raster projectExtent
#' @importFrom rgdal CRSargs
#' @importFrom gdalUtils gdalwarp
#' @importFrom raster raster
#' @importFrom rgdal GDALinfo
#' @importFrom raster intersect

EC_EnviroStack <- function(filenames, types, layernames, resamplingflag,
                           selected_layers=NULL) {
  rasters <- EC_RasterResampled (filenames, types, resamplingflag,
                                 selected_layers) # adjust to same projection, resolution and extent

  rasterstack <- raster::stack(rasters)

  names(rasterstack) <- unlist(layernames)  # assign predefined variable names
  return(rasterstack)
}

#===================================================================
### Subfunctions, listed in the order in which they appear within EC_EnviroStack

## Load rasters, determine common projection and convert categorical
## rasters to factors
EC_RasterResampled  <- function(raster_filenames, raster_types, 
                                resamplingflag, selected_layers=NULL, overwrite=TRUE) {
  rasters <- lapply(raster_filenames, EC_ReadRaster)
  
  reference <- EC_RasterRef(rasters, resamplingflag, selected_layers) # determine common raster shape
  
  rasters <- EC_RasterWarp(raster_filenames, raster_types, 
                           reference, overwrite) # adjust rasters spatially and convert categorical to factors
  
  return(rasters)
}



## Just need to be used if raster layers are not aligned up correctly
EC_RasterWarp <- function(raster_filenames, raster_types, reference, overwrite=TRUE) {
  
  rasters <- mapply(
    function(filename, type) {
      gdinfo <- rgdal::GDALinfo(filename)
      mdata <- attr(gdinfo, 'df')
      dtype <- as.character(mdata[['GDType']])
      hasNoDataValues <- mdata[['hasNoDataValue']]
      
      r <- EC_ReadRaster(filename)  # warp, crop and rescale raster file if necessary
      dir <- dirname(filename)
      temp_raster <- file.path(dir, paste0(basename(tempfile()), '.tif'))
      te <- extent(reference)
      
      # This is to fix issue with NA value being treated as value 0 if nodatavalue is not set.
      if (hasNoDataValues) {
        gdalUtils::gdalwarp(filename, temp_raster,
                            s_srs = rgdal::CRSargs(crs(r)), t_srs = rgdal::CRSargs(crs(reference)),
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
                            s_srs = rgdal::CRSargs(crs(r)), t_srs = rgdal::CRSargs(crs(reference)),
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
      if (type == "categorical") { # convert to factor if categorical
        r = as.factor(r)
      }
      return(r)
    },
    raster_filenames, raster_types)
  return(rasters)
}


## Create the same projection for all layers. If a common projection is not 
## established, it will be projected in EPSG:4326, common in EcoCommons

EC_RasterRef <- function(rasters, resamplingflag, selected_layers) {
  empty.rasters <- lapply(rasters, function(x) { raster::projectExtent(x, crs(x)) })  # create list of empty rasters
  common.crs <- crs(empty.rasters[[1]]) # choose a common.crs if all crs in rasters are the same use that one
  
  if (! do.call(raster::compareRaster, c(empty.rasters, extent=FALSE, rowcol=FALSE, 
                                 prj=TRUE, res=FALSE, orig=FALSE, rotation=FALSE, stopiffalse=FALSE))) {
    common.crs = crs("+init=epsg:4326") # common projection adopted in EcoCommons
    EC_LogWarning(sprintf("Auto projecting to common CRS: %s", common.crs))
    empty.rasters = lapply(empty.rasters, function(x) { raster::projectExtent(x, common.crs) })
  }
  
  ce <- EC_RasterExtent(empty.rasters, common.crs)
  if (! ce$equal.extents) {
    EC_LogWarning(sprintf("Auto cropping to common extent %s", EC_RasterExtentToSTR(ce$common.extent)))
  }
  
  cr <- EC_RasterResolution(empty.rasters, resamplingflag, selected_layers)
  if (! cr$is.same.res) {
    EC_LogWarning(sprintf("Auto resampling to %s resolution [%f %f]", resamplingflag,
                          cr$common.res[[1]], cr$common.res[[2]]))
  }
  
  empty.rasters <- lapply(
    empty.rasters,
    function(x) {
      extent(x) = ce$common.extent
      res(x) = cr$common.res
      return(x)
    })
  
  return(empty.rasters[[1]])  # from now on all empty.rasters should be in the same projection
}



## Apply same extension and projection for all rasters
## rasters: a vector of rasters, all rasters should have same resolution
## common.crs: crs to use to calculate intersection

EC_RasterExtent <- function(rasters, common.crs) {
  extent.list = lapply(rasters, function(r) { extent(raster::projectExtent(r, common.crs)) }) # intersect all extents
  
  common.extent = Reduce(raster::intersect, extent.list) # if all extents are the same (used to print warning)
  
  equal.extents = all(sapply(extent.list, function (x) common.extent == x))
  
  return (list(equal.extents = equal.extents, common.extent = common.extent))
}



## Extension to STR
EC_RasterExtentToSTR <- function(ext)
{
  return(sprintf("xmin=%f xmax=%f ymin=%f ymax=%f", ext@xmin, ext@xmax, ext@ymin, ext@ymax));
}



## Function to determine and resample the resolution of the user raster
## resamplingflag: a flag to determine which resampling approach to take
## selected_layers: a list of indexes to the raster layers to be considered 
## when determine the resolution to be used

EC_RasterResolution <- function(rasters, resamplingflag, selected_layers) {
  get.resolution <- function(i, rasters) {
    # Return the resolution of the raster given by the index
    return(res(rasters[[i]]))
  }
  
  if (is.null(selected_layers)) {
    resolutions = lapply(rasters, res) # Get resolutions of the raster layers
  }
  else {
    resolutions = lapply(selected_layers, get.resolution, rasters)
  }
  
  if (resamplingflag == "highest") {
    common.res = Reduce(pmin, resolutions)
  } else if (resamplingflag == "lowest") {
    common.res = Reduce(pmax, resolutions)
  }
  
  resolutions = lapply(rasters, res) # Get resolutions of all input raster layers
  is.same.res = all(sapply(resolutions, function(x) all(common.res == x)))
  return (list(common.res=common.res, is.same.res=is.same.res))
}

