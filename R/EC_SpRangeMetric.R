#' Function to compute and save species range change metric as geotif file
#'
#' @param prob_rasters 
#' @param threshold 
#' @param outfilename 
#'
#' @export
#' @importFrom 
#' 

EC_SpRangeMetric <- function(prob_rasters, threshold, outfilename) {
  # return 1 for Blank, 3 for Expansion, 0 for Contraction and 2 for No Change
  rangeChange <- overlay(as.integer(prob_rasters[[1]] >= threshold),
                         as.integer(prob_rasters[[2]] >= threshold),
                         fun=function(fp, cp) { return((2 * fp) + 1 - cp)})
  writeRaster(rangeChange, outfilename, format="GTiff",
              options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  
  # compute the area for each change category
  grid_area <- raster.from.asc(grid.area(asc.from.raster(rangeChange)))
  total_pixels = ncell(na.omit(as.data.frame(rangeChange)))
  chg_summary = as.data.frame(matrix(ncol=3, nrow=4))
  rownames(chg_summary) <- c('Contraction', 'Blank', 'No Change', 'Expansion')
  colnames(chg_summary) <- c('no_grid_cells', '%_grid_cells', 'area_km2')
  for (i in c(0,1,2,3)) {
    no_pixels = length(rangeChange[rangeChange == i])
    chg_summary[i+1,] <- c(
      no_pixels,
      (no_pixels*100.0)/total_pixels,
      sum(grid_area[rangeChange == i])/1000000.0
    )
  }
  
  # write it to a file.
  outfilename2 = outfilename
  ext = file_ext(outfilename)
  if (!is.null(ext)) {
    pattern = paste0('\\.', ext, '$')
    outfilename2 <- sub(pattern, '', outfilename)
  }
  outfilename2 = paste(outfilename2, 'csv', sep=".")
  write.csv(chg_summary, file=outfilename2, row.names=TRUE)
}