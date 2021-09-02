#' Generate a map that converts the continuous prediction into a presence/absence
#' map based on 'equal sensitivity specificity' threshold. 
#' presence-1 (red) >= threshold, absence-0 (blue) < threshold
#'
#' @param proj_raster 
#' @param threshold 
#'
#' @export EC_threshold_map
#' 

EC_threshold_map <- function(proj_raster, threshold) {
  writeLines('Generating presence/absence map ...')
  present_map <- calc(proj_raster, fun=function(x) { return(as.integer(x >= threshold))})
  
  threshold_map = file.path(EC.env$outputdir, paste0("equal.sens.spec_",
                                                     occur.species, "_maxent_cloglog.tif"))
  writeRaster(present_map, threshold_map, format="GTiff",
              options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  
  threshold_map = file.path(EC.env$outputdir, paste0("equal.sens.spec_",
                                                     occur.species, "_maxent_cloglog.png"))
  png(filename=threshold_map, width=800, height=800)
  plot(present_map, xlab='longitude', ylab='latitude')
  dev.off()
}