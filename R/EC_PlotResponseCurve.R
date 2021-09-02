#' Plot a response curve to show results - MaxEnt
#'
#' @param model 
#' @param predvars 
#' @param species 
#' @param oformat 
#'
#' @export EC_PlotResponseCurve

EC_PlotResponseCurve <- function(model, predvars, species, oformat) {  # overall response curves
  png(filename <- file.path(EC.env$outputdir, paste0("response_curves_",
                                                     species, "_maxent_", oformat, ".png")),
      width= 800, height= 800)
  response(model, at=mean, range='p', expand= 5) 
  dev.off()
  
  for (varname in predvars) {  # save individual response curve
    png(filename=file.path(EC.env$outputdir, paste0("response_curve_",
                                                    species, "_maxent_", oformat, "_", varname, ".png")),
        width= 800, height= 800)
    response(model, var=varname, at=mean, 
             fun=function(x, y, ...) predict(x, y, args=c(outputformat = oformat)), 
             range='p', expand= 5)
    dev.off()    
  }
}