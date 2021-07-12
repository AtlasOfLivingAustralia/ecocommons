#' Specify Grid Information from Geographic (lat lon) projections
#'
#' @param lats 
#' @param cellsize 
#' @param r 
#'

#' @export
#' @importFrom 
#' 

EC_GridInfo <- function(lats,cellsize, r=6378137) {
  r2 = r^2 #radius of earth

  if (length(cellsize)==1) cellsize = rep(cellsize,2) #ensure cellsize is defined for both lat & lon
  out = data.frame(lat=lats) #setup the output dataframe
  toplats = lats+(0.5 * cellsize[1]); bottomlats = lats-(0.5 * cellsize[1]) #define the top and bottom lats
  check = range(c(toplats,bottomlats), na.rm=TRUE); 
  if (-90.0001 > check[1] | 90.0001 < check[2]) 
    stop('latitudes must be between -90 & 90 inclusively')
  out$top <- EC_Distance(toplats,rep(0,length(lats)),toplats,rep(cellsize[2],length(lats)))$distance
  out$bottom <- EC_Distance(bottomlats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
  out$side <- EC_Distance(toplats,rep(0,length(lats)),bottomlats,rep(0,length(lats)))$distance
  out$diagnal <- EC_Distance(toplats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
  #calculate area of a spherical triangle using spherical excess associated by knowing distances
  #tan(E/4) = sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2))
  #where a, b, c = sides of spherical triangle
  #s = (a + b + c)/2
  excess = function(lam1,lam2,beta1,beta2){ #calculate excess... inputs are in radians
    haversine = function(y) { (1-cos(y))/2 }
    cosB1 = cos(beta1); cosB2 = cos(beta2)
    hav1 = haversine(beta2 - beta1) + cosB1 * cosB2 * haversine(lam2-lam1)
    aa = 2 * asin(sqrt(hav1)); bb = 0.5*pi - beta2; cc = 0.5*pi - beta1
    ss = 0.5 * (aa + bb + cc)
    tt = tan(ss/2) * tan((ss - aa)/2)* tan((ss - bb)/2) * tan((ss - cc)/2)
    return(abs(4 * atan (sqrt(abs(tt)))))
  }
  if (any(bottomlats== -90)) { pos = which(bottomlats== -90); 
  bottomlats[pos] = -bottomlats[pos]; toplats[pos] = -toplats[pos]} #ensure no -90 bottom lats
  out$area = excess(lam1=0,lam2=cellsize[2] *pi/180, toplats * pi/180, toplats * pi/180)
  out$area = abs(out$area - excess(lam1=0, lam2=cellsize[2] * pi/180, bottomlats * pi/180,
                                   bottomlats * pi/180))*r2
  return(out)
}