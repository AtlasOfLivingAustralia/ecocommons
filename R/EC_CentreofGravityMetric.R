#' Function to compute and save Centre of Gravity as CSV file
#'
#' @param projfiles 
#' @param outfilename 
#'
#' @export 
#' @importFrom raster raster
#' @importFrom SDMTools COGravity
#' 
#' 

EC_CentreofGravityMetric <- function(projfiles, outfilename) {
  future_proj <- raster::raster(projfiles[[1]])
  current_proj <- raster::raster(projfiles[[2]])
  future_cog <- SDMTools::COGravity(future_proj)
  current_cog <- SDMTools::COGravity(current_proj)
  
  # Do not generate CoG if it has NaN value
  if (is.nan(current_cog['COGy']) || is.nan(current_cog['COGx'])
      || is.nan(future_cog['COGy']) || is.nan(future_cog['COGx'])) {
    return()
  }
  
  results <- as.data.frame(matrix(ncol=5, nrow=3))
  rownames(results) <- c('Centre_of_Range', 'Minimum', 'Maximum')
  colnames(results) <- c('current_latitude', 'current_longitude',
                         'future_latitude', 'future_longitude', 'change_in_m')
  results[1,] = EC_Distance(current_cog['COGy'], current_cog['COGx'],
                            future_cog['COGy'], future_cog['COGx'])
  
  results[2,] = EC_Distance(min(coordinates(current_proj)[,2]),
                              min(coordinates(current_proj)[,1]),
                              min(coordinates(future_proj)[,2]),
                              min(coordinates(future_proj)[,1])
  )
 
  results[3,] = EC_Distance(max(coordinates(current_proj)[,2]),
                              max(coordinates(current_proj)[,1]),
                              max(coordinates(future_proj)[,2]),
                              max(coordinates(future_proj)[,1])
  )
  write.csv(results, file=outfilename)
}