# function to compute and save Centre of Gravity as csv file

EC_CentreofGravityMetric <- function(projfiles, outfilename) {
  future_proj = raster(projfiles[[1]])
  current_proj = raster(projfiles[[2]])
  future_cog = COGravity(future_proj)
  current_cog = COGravity(current_proj)
  
  # Do not generate CoG if it has NaN value
  if (is.nan(current_cog['COGy']) || is.nan(current_cog['COGx'])
      || is.nan(future_cog['COGy']) || is.nan(future_cog['COGx'])) {
    return()
  }
  
  results = as.data.frame(matrix(ncol=5, nrow=3))
  rownames(results) = c('Centre_of_Range', 'Minimum', 'Maximum')
  colnames(results) = c('current_latitude', 'current_longitude', 'future_latitude', 'future_longitude', 'change_in_m')
  results[1,] = BccvlDistance(current_cog['COGy'], current_cog['COGx'], future_cog['COGy'], future_cog['COGx'])
  results[2,] = BccvlDistance(min(coordinates(current_proj)[,2]),
                              min(coordinates(current_proj)[,1]),
                              min(coordinates(future_proj)[,2]),
                              min(coordinates(future_proj)[,1])
  )
  results[3,] = BccvlDistance(max(coordinates(current_proj)[,2]),
                              max(coordinates(current_proj)[,1]),
                              max(coordinates(future_proj)[,2]),
                              max(coordinates(future_proj)[,1])
  )
  write.csv(results, file=outfilename)
}