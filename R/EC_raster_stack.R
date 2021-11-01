#' Stack raster input files from user environment
#' 
#' This function adjust raster layers to same projection, resolution and extent,
#' assigning a name for the variables
#'
#' @param filenames 
#' @param types 
#' @param layernames 
#' @param resamplingflag 
#' @param selected_layers 
#' 
#' @importFrom gdalUtils gdalwarp
#' @importFrom raster compareRaster
#' @importFrom raster extent
#' @importFrom raster crs
#' @importFrom raster intersect
#' @importFrom raster projectExtent
#' @importFrom raster raster
#' @importFrom raster res
#' @importFrom raster stack
#' @importFrom rgdal CRSargs
#' @importFrom rgdal GDALinfo
#' 
#' @export EC_raster_stack


EC_raster_stack <- function (filenames,  # from predictor_info 
                             types,  # from predictor_info  
                             layernames,  # from predictor_info 
                             resamplingflag,  # from constraint_info 
                             selected_layers=NULL) {
  
  rasters <- EC_raster_resampled (filenames, types, resamplingflag,
                                 selected_layers) # adjust to same projection, resolution and extent

  rasterstack <- raster::stack(rasters)

  names(rasterstack) <- unlist(layernames)  # assign predefined variable names
  
  return(rasterstack)
}


#_____________________________________________________________________________

# Subfunctions, listed in the order in which they appear within EC_raster_stack

EC_raster_resampled  <- function (filenames,  # from predictor_info 
                                  types,  # from predictor_info  
                                  resamplingflag, # from constraint_info
                                  selected_layers=NULL,
                                  overwrite=TRUE) {
  # Load rasters, determine common projection and convert categorical to factors
  rasters <- lapply(filenames, EC_read_raster)
  
  reference <- EC_raster_ref (rasters, resamplingflag, selected_layers) # determine common raster shape
  
  rasters <- EC_raster_warp (filenames, types, 
                             reference, overwrite) # adjust rasters spatially and convert categorical to factors
  
  return(rasters)
}


#_____________________________________________________________________________


EC_raster_ref <- function (rasters,  # from EC_raster_resampled
                           resamplingflag,
                           selected_layers) {
  
  # Create the same projection for all layers. If a common projection is not 
  # established, it will be projected in WGS84 (EPSG:4326), common in EcoCommons
  empty.rasters <- lapply(rasters, function(x) { 
    raster::projectExtent(x, crs(x)) 
    })  # create list of empty rasters
  common.crs <- raster::crs(empty.rasters[[1]])  # choose a common.crs; if all crs in rasters are the same use that one
  
  if (! do.call(raster::compareRaster, c( empty.rasters,
                                          extent=FALSE, rowcol=FALSE, 
                                          prj=TRUE, res=FALSE, orig=FALSE,
                                          rotation=FALSE,
                                          stopiffalse=FALSE))) {
    
    common.crs = raster::crs("+init=epsg:4326") # common projection adopted in EcoCommons
    EC_log_warning(sprintf("Auto projecting to common CRS: %s", common.crs))
    empty.rasters = lapply(empty.rasters, function(x) { 
      raster::projectExtent(x, common.crs)
      })
  }
  
  ce <- EC_raster_extent(empty.rasters, common.crs)
  if (! ce$equal.extents) {
    EC_log_warning(sprintf("Auto cropping to common extent %s", 
                           EC_raster_to_STR(ce$common.extent)))
  }
  
  cr <- EC_raster_resolution (empty.rasters, resamplingflag, selected_layers)
  if (! cr$is.same.res) {
    EC_log_warning(sprintf("Auto resampling to %s resolution [%f %f]", resamplingflag,
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


#_____________________________________________________________________________


EC_raster_extent <- function(rasters,  # a vector of rasters, all rasters should have same resolution
                             common.crs) {  # crs to use to calculate intersection
  ## Apply same extension and projection for all rasters
  extent.list = lapply(rasters, function(r)
    { raster::extent(raster::projectExtent(r, common.crs)) }) # intersect all extents
  
  common.extent = Reduce(raster::intersect, extent.list) # if all extents are the same (used to print warning)
  
  equal.extents = all(sapply(extent.list, function (x) common.extent == x))
  
  return (list(equal.extents = equal.extents,
               common.extent = common.extent))
}


#_____________________________________________________________________________


EC_raster_to_STR <- function(ext) {
  ## Extension to STR
  return(sprintf("xmin=%f xmax=%f ymin=%f ymax=%f", ext@xmin, 
                 ext@xmax, ext@ymin, ext@ymax));
}


#_____________________________________________________________________________


EC_raster_resolution <- function (rasters,
                                  resamplingflag,  # a flag with which resampling approach to take
                                  selected_layers,
                                  reference) {  # list of indexes to the raster layers to be considered 
  # Determine and resample the resolution of the user raster  
  get.resolution <- function (i,
                              rasters) {
    # Return the resolution of the raster given by the index
    return(res(rasters[[i]]))
  }
  
  if (is.null(selected_layers)) {
    resolutions = lapply(rasters, raster::res) # Get resolutions of the raster layers
  } else {
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


#_____________________________________________________________________________


EC_raster_warp <- function (filenames,  # from predictor_info 
                            types,  # from predictor_info  
                            reference,
                            overwrite=TRUE) {
  # Just need to be used if raster layers are not aligned up correctly
  rasters <- mapply(
    function(filename, type) {
      gdinfo <- rgdal::GDALinfo(filename)
      mdata <- attr(gdinfo, 'df')
      dtype <- as.character(mdata[['GDType']])
      hasNoDataValues <- mdata[['hasNoDataValue']]
      
      r <- EC_read_raster(filename)  # warp, crop and rescale raster file if necessary
      dir <- dirname(filename)
      temp_raster <- file.path(dir, paste0(basename(tempfile()), '.tif'))
      te <- extent(reference)
      
      # Fix issue with NA value being treated as value 0 if nodatavalue is not set.
      if (hasNoDataValues) {
        gdalUtils::gdalwarp(filename, temp_raster,
                            s_srs = rgdal::CRSargs(crs(r)),
                            t_srs = rgdal::CRSargs(crs(reference)),
                            te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                            ts=c(ncol(reference), nrow(reference)),
                            # tr=c(...), ... either this or ts
                            r="near",
                            of="GTiff",
                            dstnodata=mdata[['NoDataValue']],
                            co=c("TILED=YES", "COMPRESS=LZW"))
      } else {
        gdalUtils::gdalwarp(filename, temp_raster,
                            s_srs = rgdal::CRSargs(crs(r)),
                            t_srs = rgdal::CRSargs(crs(reference)),
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
    
    filenames, types)
  
  return(rasters)
}