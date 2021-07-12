#' Function created to get family from a string (character) in a safe way
#' works for all variants of the R family object (e.g. see ?family)
#'
#' @param s 
#'
#' @export EC_FamilyFromString
#' @ImportFrom
#' 

EC_FamilyFromString <- function(s) {
  s = gsub(pattern="\"|| ", replacement="", s) # strip quotes and spaces
  f = gsub(pattern="\\(.*\\)", replacement="", s) # the name of the function
  
  allowable= c("binomial",
               "gaussian",
               "Gamma",
               "inverse.gaussian",
               "poisson",
               "quasi",
               "quasibinomial",
               "quasipoisson")
  
  if (! f %in% allowable )
  {
    stop(sprintf("unsupported function %s", f))
  }
  
  fargs=gsub(pattern=".*\\(||\\)",
             replacement = "",
             sub(pattern = f,
                 replacement = "",
                 s)) #get the args inside the parentheses
  args=list()
  
  if (fargs != "")
  {
    l=strsplit(fargs, ",")[[1]]
    for( i in 1:length(l) )
    {
      ll=strsplit(l[i],"=")[[1]]
      if (length(ll) == 2)
      {
        args[ll[1]] = ll[2]
      }
      else
      {
        stop(sprintf("unhandled result when splitting %s", l[i]))
      }
    }
  }
  return (do.call(what=f, args=args))
}
