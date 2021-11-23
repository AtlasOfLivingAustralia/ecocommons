#' Function to geographically constrained the Species Distribution Modelling.
#' The constraint is the intersection between the occurrence's convex-hull 
#' polygon and the constraint. Otherwise, actual constraint is the convex-hull 
#' polygon
#' A buffer of width 1-resolution cell is applied to avoid missing environment
#' values along the boundary of the polygon
#'
#' @param rasterstack 
#' @param occur 
#' @param absen 
#' @param constraints 
#' @param generateCHull 
#' 
#' @importFrom raster compareCRS
#' @importFrom raster crop
#' @importFrom raster mask
#' @importFrom raster stack
#' @importFrom rgeos gBuffer
#' @importFrom rgeos gConvexHull
#' @importFrom rgdal readOGR
#' @importFrom rgdal writeOGR
#' @importFrom rjson fromJSON
#' @importFrom sp SpatialPolygons
#' @importFrom sp SpatialPoints
#' @importFrom sp spTransform
#' 
#' @export EC_SDM_geoconstrained

EC_SDM_geoconstrained <- function(current_climate, 
                                 occur, 
                                 absen, 
                                 constraints, 
                                 generateCHull) {

  if (is.null(constraints)) {
    parsed_geojson <- sp::SpatialPolygons(list(Polygons(list(Polygon(rbind(c(1,1)))), ID=1)),
                                      proj4string=crs('+init=epsg:4326'))
  } else {
    parsed_geojson <- rgdal::readOGR(dsn = paste0(constraints),
                                     layer = "OGRGeoJSON", verbose = FALSE, 
                                     p4s = '+proj=longlat +datum=WGS84 +no_defs')
  }
  
  if (!compareCRS(current_climate, parsed_geojson, verbatim=TRUE)) {
    crs(current_climate) <-crs(parsed_geojson)  # if CRS is different, re-project
  }
  
  occurSP <- sp::SpatialPoints(occur,  #transform your occurrence data in a Polygon (spatial points)
                                 proj4string=crs('+init=epsg:4326'))  # put in the same projection as the layers
    
 # constraint is the intersection between the occurrence's convex-hull polygon and the constraint
  region_offset = 0
    
 # otherwise, actual constraint is the convex-hull polygon
  if (generateCHull) {
    # get the offset 
    constraint_json <- rjson::fromJSON(constraints)
    region_offset <- constraint_json$properties$region_offset
    if (is.null(region_offset) || is.na(region_offset) || region_offset == '') {
      region_offset = 0
      } else {
      region_offset <- as.double(region_offset)
      region_offset <- ifelse(!is.na(region_offset) && is.numeric(region_offset),
                                region_offset/111.0, 0) # convert from km to degree
      }
      
    chullPolygon <- rgeos::gConvexHull(occurSP)  # generate a convex hull polygon with species distribution data
      
    if (!is.null(constraints)) {
        parsed_geojson <- intersect(parsed_geojson, chullPolygon)  # intersect the constaint region with the convex hull polygon 
        } else {
          parsed_geojson <- chullPolygon
        }
    }
  
  parsed_geojson <- rgeos::gBuffer(parsed_geojson, byid=TRUE, width=max(region_offset,
                                             max(res(current_climate@layers[[1]])))) # create buffer of width 1-resolution cell
    
  if (generateCHull) {
    parsed_geojson_df <- data.frame(ID=1:length(parsed_geojson))
    parsed_geojson_df <- SpatialPolygonsDataFrame(parsed_geojson, parsed_geojson_df) 
    filename <- file.path(EC.env$outputdir, 'modelling_region.json')
    rgdal::writeOGR(parsed_geojson_df, dsn = filename, layer = 'parsed_geojson_df',
                      driver='GeoJSON')  # save the convex-hull generated as geojson
    
    # read the json to save the projection (crs) and properties from original to the new json  
    gjson <- rjson::fromJSON(file=filename)
    origjson <- rjson::fromJSON(constraints)
    
    if (! 'crs' %in% names(gjson)) {
      gjson = append(gjson, list(crs=origjson$crs))   
    }
    
    if (! 'properties' %in% names(gjson)) {
      gjson = append(gjson, list(properties=origjson$properties))
    }
    
    write(rjson::toJSON(gjson), filename)

    origjson = NULL
    gjson = NULL
  }
    
  occurSPconstrained <- occurSP[!is.na(over(occurSP, parsed_geojson))]  # remove NAs and points that don't overlap with constrained area
  occurconstrained <- as.data.frame(occurSPconstrained)
    
  # constraint the true absence points, if available
  absenconstrained = NULL
  
  if (!is.null(absen) && nrow(absen) > 0) {
    absenSP <- EC_data_projection(absen, current_climate)
    absenSPconstrained <- absenSP[!is.na(over(absenSP, parsed_geojson))]  # remove NAs and points that don't overlap with constrained area
      
    absenSPconstrained <- spTransform(absenSPconstrained, crs('+init=epsg:4326')) # project it back to epsg:4326 for saving as a csv file
    absenconstrained <- as.data.frame(absenSPconstrained)
    names(absenconstrained) <- c("lon", "lat")
    write.csv(absen, file=absenFilename, row.names=FALSE)
  }
  
  
  # crop and mask the rasterstack (and make sure it is a RasterStack)
  cropped_raster <- raster::crop(current_climate,
                                 extent(parsed_geojson))  # crop to constraint region before masking
  
  envraster_filename <- paste(EC.env$outputdir, 
                              basename(tempfile(fileext = ".grd")), sep = "/")
  
  geoconstrained <- raster::mask(cropped_raster, parsed_geojson,
                                 filename = envraster_filename)
  
  geoconstrained <- stack(geoconstrained)
  
  EC_raster_remove(stack(cropped_raster))
  
  # return the masked raster stack and constrained occurrence points
  
  mylist <- list("masked_raster" = geoconstrained,
                 "occur" = occurconstrained,
                 "absen" = absenconstrained)
  
  return(mylist)
  
}




#_____________________________________________________________________________
## Subfunction, as appear within EC_SDM_geoconstrained()


EC_data_projection <- function(data,
                               climate_data) {
  # Project species distribution data and climate data in the same crs (EPSG:4326)
  sp <- SpatialPoints(data)
  if (is.na(crs(sp))) {
    crs(sp) <- '+init=epsg:4326'
  }
  
  if (!raster::compareCRS(sp, climate_data, verbatim=TRUE)) {
    sp <- spTransform(sp, crs(climate_data))
  }
  return(sp)
}
