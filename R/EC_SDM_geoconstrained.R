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
  
  if (is.null(constraints) & !generateCHull) {
    return(list("raster_obj" = current_climate, "occur" = occur,
                "absen" = absen))
  }
  
  geojson_crs <- CRS("+init=epsg:4326") # create a geojson for convex-hull polygon if no geojson
  if (is.null(constraints)) {
    parsed_geojson <- sp::SpatialPolygons(list(Polygons(list(Polygon(rbind(c(1,1)))), ID=1)),
                                      proj4string=crs('+init=epsg:4326'))
  } else {
    parsed_geojson <- rgdal::readOGR (dsn = paste0(constraints),
                                     layer = "OGRGeoJSON", verbose = FALSE)
    geojson_crs <- crs(parsed_geojson)
  }
  
  if (!compareCRS(current_climate, parsed_geojson, verbatim=TRUE)) {
    crs(current_climate) <-crs(parsed_geojson)  # if CRS is different, re-project
  }
  
  if (!is.null(occur)) {
    
    occurSP <- sp::SpatialPoints(occur)  # constrain the occurrence points
    if (is.na(crs(occurSP))) {
      crs(occurSP) <- '+init=epsg:4326'
    }
    
    if (!compareCRS(occurSP, parsed_geojson, verbatim=TRUE)) {
      occurSP <- sp::spTransform(occurSP, crs(parsed_geojson))
    }
    
    region_offset = 0
    if (generateCHull) {
      # get the offset 
      constraint_json <- rjson::fromJSON(constraints)
      region_offset <- constraint_json$properties$region_offset
      if (is.null(region_offset) || is.na(region_offset) || region_offset == '') {
        region_offset = 0
      }
      else {
        region_offset <- as.double(region_offset)
        region_offset <- ifelse(!is.na(region_offset) && is.numeric(region_offset),
                                region_offset/111.0, 0) # convert from km to degree
      }
      
      chcoords <- occurSP@coords[chull(occurSP@coords[,1:2]),]
      chullPolygon <- sp::SpatialPolygons(list(Polygons(list(Polygon(chcoords[,1:2],
                                                                         hole=FALSE)),
                                                    ID=1)), proj4string=crs(parsed_geojson))
      if (!is.null(constraints)) {
        parsed_geojson <- intersect(parsed_geojson, chullPolygon)
      }
      else {
        parsed_geojson <- chullPolygon
      }
    }
    
    parsed_geojson <- rgeos::gBuffer(parsed_geojson, byid=TRUE, width=max(region_offset,
                                               max(res(current_climate@layers[[1]])))) #create buffer of width 1-resolution cell
    
    if (generateCHull) {
      filename <- file.path(EC.env$outputdir, 'modelling_region.json')
      transformed_geojson <- sp::spTransform(parsed_geojson, geojson_crs) 
      rgdal::writeOGR(transformed_geojson, filename, 'OGRGeoJSON', driver='GeoJSON')  # save the convex-hull generated as geojson.
      transformed_geojson = NULL
      
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
    
    occurSPconstrained <- occurSP[!is.na(over(occurSP, parsed_geojson))]
    occurconstrained <- as.data.frame(occurSPconstrained)
    
    absenconstrained = NULL

    if (!is.null(absen) && nrow(absen) > 0) {
      absenSP <- EC_data_projection(absen, current_climate)
      absenSPconstrained <- absenSP[!is.na(over(absenSP, parsed_geojson))]
      
      absenSPconstrained <- spTransform(absenSPconstrained, CRS('+init=epsg:4326')) # project it back to epsg:4326 for saving as a csv file
      absenconstrained <- as.data.frame(absenSPconstrained)
      names(absenconstrained) <- c("lon", "lat")
      write.csv(absen, file=absenFilename, row.names=FALSE)
    } else {
    occurconstrained = NULL
    absenconstrained = NULL
  }
  
  geoconstrained <- raster::stack(current_climate, EC_mask, parsed_geojson)  # mask the rasterstack (current_climate)
  
  mylist <- list("raster_obj" = geoconstrained, "occur" = occurconstrained,
                 "absen" = absenconstrained)
  return(mylist)
}




#_____________________________________________________________________________
## Subfunction, as appear within EC_SDM_geoconstrained()


EC_data_projection <- function(data,
                               climate.data) {
  # Project species distribution data and climate data in the same crs (EPSG:4326)
  sp <- SpatialPoints(data)
  if (is.na(crs(sp))) {
    crs(sp) <- '+init=epsg:4326'
  }
  
  if (!raster::compareCRS(sp, climate.data, verbatim=TRUE)) {
    sp <- spTransform(sp, crs(climate.data))
  }
  return(sp)
  }
}


#_____________________________________________________________________________


EC_mask <- function(current_climate,
                    parsed_geojson) {
  # Crop the raster to the extent of the constraint region and mask it
  cropped_raster <- raster::crop (current_climate,
                                  extent(parsed_geojson))  # crop to constraint region before masking
  
  envraster_filename <- paste(EC.env$workdir, 
                              basename(tempfile(fileext = ".grd")), sep="/")
  
  masked_raster <- raster::mask (cropped_raster, parsed_geojson,
                                 filename = envraster_filename)
  
  if (is.factor(masked_raster)) {
    masked_raster = as.factor(masked_raster)
  }
  
  EC_raster_remove(stack(cropped_raster))  # remove cropped raster and associated raster files (i.e. grd and gri)
  
  return(masked_raster)
}
