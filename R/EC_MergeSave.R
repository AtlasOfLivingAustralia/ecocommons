# This is a common function to merge and save the new generated data
# in a csv file

EC_MergeSave <- function(env, csvdata, spname, ofname) {
  data <- cbind(csvdata, species=spname, extract(env, csvdata[c('lon','lat')]))
  
  EC_WriteCSV(data, ofname, rownames=FALSE)
}