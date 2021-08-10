# Function to geographically constrained the species distribution modelling
#' The constraint is the intersection between the occurrence's convex-hull 
#' polygon and the constraint. Otherwise, actual constraint is the convex-hull 
#' polygon
#' A buffer of width 1-resolution cell is applied to avoid missing environment
#' values along the boundary of the polygon
#'
#' @param rasterstack 
#' @param occur 
#' @param absen 
#' @param raw_geojson 
#' @param generateCHull 
#'
#' @export EC_SDMGeoConstrained
#' @importFrom rgdal readOGR
#' @importFrom rgeos gBuffer
#' @importFrom rjson fromJSON
#' @importFrom sp SpatialPolygons
#' @importFrom sp SpatialPoints
#' @importFrom raster compareCRS
#' @importFrom raster crop
#' @importFrom raster mask
#'            
#'

EC_SDMGeoConstrained <- function(rasterstack, occur, absen, raw_geojson, generateCHull) {
  
  if (is.null(raw_geojson) & !generateCHull) {
    return(list("raster" = rasterstack, "occur" = occur, "absen" = absen))
  }
  
  geojson_crs <- CRS("+init=epsg:3857") # create a geojson for convex-hull polygon if no geojson
  if (is.null(raw_geojson)) {
    parsed_geojson <- sp::SpatialPolygons(list(Polygons(list(Polygon(rbind(c(1,1)))), ID=1)),
                                      proj4string=crs(rasterstack))
  } else {
    parsed_geojson <- rgdal::readOGR(dsn = raw_geojson, layer = "OGRGeoJSON", verbose = FALSE)
    geojson_crs <- crs(parsed_geojson)
  }
  
  if (!compareCRS(rasterstack, parsed_geojson, verbatim=TRUE)) {
    parsed_geojson <- spTransform(parsed_geojson, crs(rasterstack))  # if CRS is different, reproject geojson to rasterstack
  }
  
  if (!is.null(occur)) {
    
    occurSP <- sp::SpatialPoints(occur)  # constrain the occurrence points
    if (is.na(crs(occurSP))) {
      crs(occurSP) <- '+init=epsg:4326'
    }
    
    if (!compareCRS(occurSP, parsed_geojson, verbatim=TRUE)) {
      occurSP <- spTransform(occurSP, crs(parsed_geojson))
    }
    
    region_offset = 0
    if (generateCHull) {
      # get the offset 
      constraint_json <- rjson::fromJSON(raw_geojson)
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
      chullPolygon <- raster::SpatialPolygons(list(Polygons(list(Polygon(chcoords[,1:2],
                                                                         hole=FALSE)),
                                                    ID=1)), proj4string=crs(parsed_geojson))
      if (!is.null(raw_geojson)) {
        parsed_geojson <- intersect(parsed_geojson, chullPolygon)
      }
      else {
        parsed_geojson <- chullPolygon
      }
    }
    
    parsed_geojson <- rgeos::gBuffer(parsed_geojson, byid=TRUE, width=max(region_offset,
                                                                   max(res(rasterstack@layers[[1]])))) #create buffer of width 1-resolution cell
    
    if (generateCHull) {
      filename <- file.path(EC.env$outputdir, 'modelling_region.json')
      transformed_geojson <- spTransform(parsed_geojson, geojson_crs) 
      writeOGR(transformed_geojson, filename, 'OGRGeoJSON', driver='GeoJSON')  # save the convex-hull generated as geojson.
      transformed_geojson = NULL
      
      gjson <- rjson::fromJSON(file=filename)
      origjson <- rjson::fromJSON(raw_geojson)
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
      absenSP <- EC_SpDataProjection(absen, rasterstack)
      absenSPconstrained <- absenSP[!is.na(over(absenSP, parsed_geojson))]
      
      absenSPconstrained <- spTransform(absenSPconstrained, CRS('+init=epsg:4326')) # project it back to epsg:4326 for saving as a csv file
      absenconstrained <- as.data.frame(absenSPconstrained)
      #names(absenconstrained) <- c("lon", "lat")
      #write.csv(absen, file=absenFilename, row.names=FALSE)
    }
  }
  else {
    occurconstrained = NULL
    absenconstrained = NULL
  }
  
  geoconstrained <- stack(lapply(as.list(rasterstack), EC_Mask, parsed_geojson))  # mask the rasterstack
  
  mylist <- list("raster" = geoconstrained, "occur" = occurconstrained, "absen" = absenconstrained)
  return(mylist)
}




#________________________________________________________________
## Subfunctions, listed in the order in which they appear within EC_SDMGeoConstrained()


## Project species distribution data and climate data in the same crs (EPSG:4326)
EC_SpDataProjection <- function(data, climate.data) {
  sp <- SpatialPoints(data)
  if (is.na(crs(sp))) {
    crs(sp) <- '+init=epsg:4326'
  }
  
  if (!raster::compareCRS(sp, climate.data, verbatim=TRUE)) {
    sp <- spTransform(sp, crs(climate.data))
  }
  return(sp)
}




## Function to crop the raster to the extent of the constraint region and 
## mask it

EC_Mask <- function(raster, parsed_geojson) {
  
  cropped_raster <- raster::crop(raster, extent(parsed_geojson))  # crop to constraint region before masking
  
  envraster_filename <- paste(EC.env$workdir, basename(tempfile(fileext = ".grd")), sep="/")
  masked_raster <- raster::mask(cropped_raster, parsed_geojson, filename = envraster_filename)
  
  if (is.factor(masked_raster)) {
    masked_raster = as.factor(masked_raster)
  }
  
  EC_RevRasterObject(stack(cropped_raster))  # remove cropped raster and associated raster files (i.e. grd and gri)
  
  return(masked_raster)
}
