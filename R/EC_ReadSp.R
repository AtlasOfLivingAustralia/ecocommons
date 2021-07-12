# Read species presence/absence data
#
# To run Species Distribution Models (SDMs) we first need to set up the
# necessary dataset. This script aims to check if the species geographic
# data is on the right format to be read. It means that latitude and longitude
# should be numeric, and month and year information shoud be kept.

EC_ReadSp <- function(filename, month_filter=NULL) {
  if (!is.null(filename)) {  # return NULL if filename is not given
    csvfile <- read.csv(filename, colClasses=c("lon"="numeric", "lat"="numeric"))

    col_kept <- c("lon","lat") # keep only lon and lat columns; for MM include month column
    if ('year' %in% colnames(csvfile)) {
      col_kept = append(col_kept, 'year') # keep year column if exists
    }
    if (is.null(month_filter)) {
      csvfile = csvfile[col_kept]
      return(csvfile)
    }
    else {
      col_kept = append(col_kept, 'month')
      csvfile = csvfile[col_kept]
      return(subset(csvfile, month %in% unlist(month_filter)))
    }
  }
}
