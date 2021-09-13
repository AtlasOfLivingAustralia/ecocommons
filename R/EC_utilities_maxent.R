#' Subfunctions ONLY to run the MaxEnt algorithm. It builds evaluation parameters,
#' set data on biomod2 structure and save the result
#'
#' @importFrom dismo randomPoints
#' @importFrom raster extract
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom sp SpatialPoints
#'
#' Not exported

EC_create_background <- function(rastermask,
                                 bgsize) {
  # Generate background data on the region of interest
  background <- dismo::randomPoints(rastermask, n=bgsize,
                                    extf=1.0, excludep=FALSE)

  colnames(background) <- c('lon', 'lat')

  return(as.data.frame(background))
}


#_____________________________________________________________________________


EC_create_biased_background <- function(bg,
                                        biasfile,
                                        bgsize) {
  # If chosen to use a bias file, apply it in the occurences to the background
  # points
  bias <- raster::raster(biasfile)

  bg_sp_class<-sp::SpatialPoints(coords = cbind(bg$lon,bg$lat),
                                 proj4string = rgdal::CRS(as.character(NA)))   # convert background points to spatial points formal class

  weight <- raster::extract(bias, bg_sp_class)

  bg_sp_class <- bg_sp_class[!is.na(weight), ]

  weight <- weight[!is.na(weight)]  # random sample of the candidate background points

  bg_biased <- bg_sp_class[sample(length(bg_sp_class), size=bgsize,
                                  replace = TRUE, prob=weight), ]

  colnames(bg_biased@coords) <- c("lon", "lat")

  bg_biased <- bg_biased@coords

  return(bg_biased)
}


#_____________________________________________________________________________



EC_plot_contin_OCC <- function (dataToExplore,
                                listVars,
                                doCorrelation = TRUE,
                                outerTitle = "",
                                fnamePrefix = "occurrence_predictors_",
                                outputdir = "output")  {
  # Load and export plot function for continuous predictor variables for
  # the occurrence data

  nVars = length(listVars)

  if (nVars == 0) {
    return()
  }

  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)
  plotWidth = 800; plotHeight = 800

  # Boxplots
  png(filename=file.path(outputdir, paste0(fnamePrefix, 'boxplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    boxplot(dataToExplore[,c(nameVar)], main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

  # Histograms
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'histogram.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    hist(dataToExplore[,c(nameVar)], breaks = 10, xlab = "", main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

  # Density plots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'density.png')),
      width =plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    plot(density(dataToExplore[,c(nameVar)], na.rm = TRUE), main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

  # Pearson Correlation
  if (doCorrelation && nCols > 1) {

    png(filename = file.path(outputdir, paste0(fnamePrefix, 'correlations.png')),
        width = plotWidth, height = plotHeight)
    par(oma = c(0, 0, 1.5, 0), ps = 15)
    pairs(dataToExplore,
          lower.panel = function(x, y){
            usrSaved = par("usr"); on.exit(par(usrSaved))
            par(usr = c(0, 1, 0, 1), ps = 15)
            text(0.5, 0.5, cex = 2,
                 round(cor(x, y, use="na.or.complete", method="pearson"), digits = 2))})
    title(outerTitle, outer = TRUE, line = 0)
    dev.off()
  }
}


#_____________________________________________________________________________


EC_plot_contin_BG <- function (dataToExplore,
                               listVars,
                               doCorrelation = TRUE,
                               outerTitle = "",
                               fnamePrefix = "background_predictors_",
                               outputdir = "output")  {
  # Load and export plot function for continuous predictor variables for the
  # background data
  nVars = length(listVars)

  if (nVars == 0) {
    return()
  }

  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)
  plotWidth = 800; plotHeight = 800

  # Boxplots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'boxplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    boxplot(dataToExplore[,c(nameVar)], main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

  # Histograms
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'histogram.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    hist(dataToExplore[,c(nameVar)], breaks = 10, xlab = "", main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

  # Density plots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'density.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    plot(density(dataToExplore[,c(nameVar)], na.rm = TRUE), main = nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()
}


#_____________________________________________________________________________


EC_plot_categ_OCC <- function(dataToExplore,
                              listVars,
                              outerTitle = "",
                              fnamePrefix = "occurrence_predictors_",
                              outputdir = "output") {
  # Load and export plot function for categorical predictor variables for
  # the occurrence data

  nVars = length(listVars)

  if (nVars == 0) {
    return()
  }

  plotWidth = 800; plotHeight = 800
  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)

  # Barplots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'barplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows, nCols), oma = c(0, 0, 1.5, 0), ps = 15)
  for (nameVar in listVars) {
    barplot(table(dataToExplore[,c(nameVar)]),
            main=nameVar)
  }
  title(outerTitle, outer = TRUE, line = 0)
  dev.off()

}


#_____________________________________________________________________________


EC_plot_categ_BG <- function(dataToExplore,
                             listVars,
                             outerTitle = "",
                             fnamePrefix = "background_predictors_",
                             outputdir = "output") {
  # Load and export plot function for categorical predictor variables for
  # the background data

  nVars = length(listVars)

  if (nVars == 0) {
    return()
  }

  plotWidth = 800; plotHeight = 800
  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)

  # Barplots

  png(filename = file.path (outputdir, paste0(fnamePrefix, 'barplot.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps = 15)

  for (nameVar in listVars) {
    barplot(table(dataToExplore[,c(nameVar)]),
            main = nameVar)
  }

  title(outerTitle, outer = TRUE, line = 0)

  dev.off()
}


#_____________________________________________________________________________


EC_plot_response_curve <- function (model,
                                    predvars,
                                    species,
                                    oformat) {  # overall response curves

  # Response curve to show maxEnt results
  png(filename <- file.path(EC.env$outputdir, paste0("response_curves_",
                                                     species, "_maxent_",
                                                     oformat, ".png")),
      width= 800, height= 800)

  response(model, at = mean, range = 'p', expand = 5)
  dev.off()

  for (varname in predvars) {  # save individual response curve
    png(filename <- file.path(EC.env$outputdir, paste0("response_curve_",
                                                       species, "_maxent_",
                                                       oformat, "_",
                                                       varname, ".png")),
        width= 800, height= 800)
    response(model, var = varname, at = mean,
             fun = function(x, y, ...) predict(x, y, args= c(outputformat = oformat)),
             range='p', expand= 5)
    dev.off()
  }
}


#_____________________________________________________________________________


EC_threshold_map <- function(proj_raster,
                             threshold) {

  # Generate a map that converts the continuous prediction into a presence/absence
  # map based on 'equal sensitivity specificity' threshold.
  # presence-1 (red) >= threshold, absence-0 (blue) < threshold

  writeLines('Generating presence/absence map ...')
  present_map <- calc(proj_raster, fun=function(x) {
    return(as.integer(x >= threshold))})

  threshold_map = file.path(EC.env$outputdir,
                            paste0("equal.sens.spec_",
                                   occur.species, "_maxent_cloglog.tif"))

  raster::writeRaster (present_map, threshold_map, format="GTiff",
                       options=c("COMPRESS=LZW", "TILED=YES"),
                       overwrite=TRUE)

  threshold_map = file.path(EC.env$outputdir, paste0("equal.sens.spec_",
                                                     occur.species, "_maxent_cloglog.png"))

  png(filename=threshold_map, width=800, height=800)

  plot(present_map, xlab='longitude', ylab='latitude')

  dev.off()
}
