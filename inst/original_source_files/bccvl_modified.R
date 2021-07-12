############################################################
###        EcoCommons script to define parameters        ###
############################################################
##
## Author details: EcoCommons Platform, Contact details: emailadress@
## Copyright statement: This script is the product of BCCVL etc.
## Date : ??
## Script and data info:
## In this script you will create functions to set up parameters that
## will allow you to run multiple algorithms for species distribution
## modelling (SDMs). To start, you will need 1) install and load all
## necessary packages, 2) load a document with parameters?


# First, setup a R environment on host computer, defining a CRAN repository
# You can jump this step if running remotely at home
r <- getOption ("repos")
r["CRAN"] <- "http://cran.ms.unimelb.edu.au/"
options(repos = r)
# print warnings immediately
options(warn = 1)


# Check the required packages for distribution models
# 1. print out list of installed packages
write.table(installed.packages() [,c("Package", "Version", "Priority")],
            row.names=FALSE)

# 2. check if libraries are installed, install if necessary and then load them
necessary <- c("ggplot2","tools", "rjson", "dismo", "gbm", "rgdal", "rgeos", "pROC", "png", "gstat", "gdalUtils", "ggdendro", "raster","biomod2","rasterVis") #list the libraries needed

###REMOVED NOW JUST FOR TEST:::: "SDMTools",‘spatial.tools’ and ‘rmaxent - available for R 4.0 on github’


installed <- necessary %in% installed.packages() #check if library is installed
if (length(necessary[!installed]) >= 1) {
    install.packages(necessary[!installed], dep = T) #if library is not installed, install it
}
for (lib in necessary) {
    library(lib,character.only = T) #load the libraries
}

# Load the parameters
# The 'param.json' file used here is generated on EcoCommons python script
params <- rjson::fromJSON(file="~/Documents/BCCVL_scripts/test_script.json")
EC.params <- params$params
EC.env <- params$env
rm(params)

# Set working directory (script runner takes care of it)
setwd(EC.env$outputdir)
# Set a temporary directory to raster - we do this because raster sometimes makes
# temp files (e.g. when cropping).Might want to make this configurable - e.g. we
# might want to control maxmemory, and/or other raster options
rasterOptions(tmpdir=paste(EC.env$workdir,"raster_tmp",sep="/"))

# Use seed supplied, if any. Otherwise generate a random seed.
seed <- EC.params$random_seed
if (is.null(seed)) {
    seed = runif(1, -2^31, 2^31-1)
}
seed <- as.integer(seed)
set.seed(seed)
EC.params["random_seed"] <- seed

############################################################
# Define helper functions to print parameters
############################################################

a<-EC_RenameParameters(EC.params)

EC_RenameParameters <- function (param, value) {
  # Computes simple names for different parameters
  #
  # Args:
  #  param: list of possible parameters that can be implemented
  #  value: ?
  #
  # Return:
  #  pname: new names for the different parameters, including value
    pname <- gsub("_", " ", param)
    if (param == "prevalence") {
        pname = "weighted response weights"
    }
    else if (param == "var_import") {
        pname = "resampling"
    }
    else if (param == "nbcv") {
        pname = "NbCV"
    }
    else if (param == "n_trees") {
        pname = "trees added each cycle"
    }
    else if (param == "control_xval") {
        pname = "cross-validations"
    }
    else if (param == "control_minbucket") {
        pname = "minimum bucket"
    }
    else if (param == "control_minsplit") {
        pname = "minimum split"
    }
    else if (param == "control_cp") {
        pname = "complexity parameter"
    }
    else if (param == "control_maxdepth") {
        pname = "maximum depth"
    }
    else if (param == "irls_reg") {
        pname = "irls.reg"
    }
    else if (param == "maxit") {
        pname = "maximum iterations"
    }
    else if (param == "mgcv_tol") {
        pname = "convergence tolerance"
    }
    else if (param == "mgcv_half") {
        pname = "number of halvings"
    }
    else if (param == "n_minobsinnode") {
        pname = "Min observations in terminal node"
    }
    else if (param == "control_epsilon") {
        pname = "control: epsilon"
    }
    else if (param == "control_maxit") {
        pname = "control: maxit"
    }
    else if (param == "control_trace") {
        pname = "control: trace"
    }
    else if (param == "model") {
        pname = "Model returned"
    }
    else if (param == "x") {
        pname = "x returned"
    }
    else if (param == "y") {
        pname = "y returned"
    }
    else if (param == "qr") {
        pname = "QR returned"
    }
    else if (param == "singular_ok") {
        pname = "Singular fit ok"
    }
    else if (param == "thresh") {
        pname = "threshold"
    }
    else if (param == "maximumiterations") {
        pname = "Maximum iterations"
    }
    else if (param == "ntree") {
        pname = "number of trees"
    }
    else if (param == "mtry") {
        pname = "number of variables at each split (mtry)"
    }
    else if (param == "nodesize") {
        pname = "node size"
    }
    else if (param == "maxnodes") {
        pname = "maximum nodes"
    }
    else if (param == "pa_ratio") {
        pname = "absence-presence ratio"
    }
    return(paste(pname, " = ", value, "\n", sep="", collapse=""))
}

EC.ParameterPrint <- function(params) {
  # Assign parameters necessary to run specific algorithms
  #
  # Args:
  #  parameters: set of parameters that can be used
  #  func: functions or algorithms that the parameters are assigned to
  #
  # Return:
  #  list of parameters restricted to specific algorithms
   func <- params[["function"]]
    if (is.null(func))
        return("")
    cat("Algorithm:", func, "\n")

    pnames = c("random_seed")
    if (func == "ann") {
        pnames = c("prevalence", "var_import", "maxit", "nbcv", "rang", "random_seed")
    }
    else if (func == "brt") {
        pnames = c("tree_complexity", "learning_rate", "bag_fraction", "n_folds", "prev_stratify", "family", "n_trees", "max_trees", "tolerance_method", "tolerance_value", "random_seed")
    }
    else if (func == "cta") {
        pnames = c("prevalence", "var_import", "method", "control_xval", "control_minbucket", "control_minsplit", "control_cp", "control_maxdepth", "random_seed")
    }
    else if (func == "fda") {
        pnames = c("prevalence", "var_import", "method", "random_seed")
    }
    else if (func == "gam") {
        pnames = c("prevalence", "var_import", "interaction_level", "family", "irls_reg", "epsilon", "maxit", "mgcv_tol", "mgcv_half", "random_seed")
    }
    else if (func == "gamlss") {
        pnames = c("sigma_formula", "nu_formula", "tau_formula", "family", "weights", "contrasts", "method", "start_from", "mu_start", "sigma_start", "nu_start", "tau_start", "mu_fix", "sigma_fix", "nu_fix", "tau_fix", "control", "i_control", "other_args", "random_seed")
    }
    else if (func == "gbm") {
        pnames = c("prevalence", "var_import", "distribution", "n_trees", "interaction_depth", "n_minobsinnode", "shrinkage", "bag_fraction", "train_fraction", "cv_folds", "random_seed")
    }
    else if (func == "glm") {
        pnames = c("prevalence", "var_import", "type", "interaction_level", "test", "family", "mustart", "control_epsilon", "control_maxit", "control_trace", "random_seed")
    }
    else if (func == "lm") {
        pnames = c("subset", "weights", "na_action", "method", "model", "x", "y", "qr", "singular_ok", "contrasts", "offset", "random_seed")
    }
    else if (func == "manova") {
        pnames = c("projections_returned", "qr", "contrasts", "subset", "weights", "na_action", "random_seed")
    }
    else if (func == "mars") {
        pnames = c("prevalence", "var_import", "degree", "nk", "penalty", "thresh", "prune", "random_seed")
    }
    else if (func == "maxent") {
        pnames = c("prevalence", "var_import", "maximumiterations", "linear", "quadratic", "product", "threshold", "hinge", "lq2lqptthreshold", "lq2lqthreshold", "hingethreshold", "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "defaultprevalence", "random_seed")
    }
    else if (func == "rf") {
        pnames = c("prevalence", "var_import", "do.classif", "ntree", "mtry", "nodesize", "maxnodes", "random_seed")
    }
    else if (func == "sre") {
        pnames = c("prevalence", "var_import", "quant", "random_seed")
    }

    pnames = c(pnames, "pa_ratio", "pa_strategy", "pa_sre_quant", "pa_disk_min", "pa_disk_max")
    for (p in pnames) {
        cat(EC.RenameParameters(p, params[[p]]))
    }
    return("")
}


# Print out parameters used
EC.ParameterPrint(EC.params)

############################################################
# Define helper functions to use in EC
############################################################

## Needed for catching
EC_ErrNull <- function (e) return (NULL)

# Read species presence/absence data
# Create a function - return NULL if filename is not given
EC_SpRead <- function(filename, month_filter=NULL) {
    if (!is.null(filename)) {
        # We might loose precision of lon/lat when onverting to double,
        # however, given the nature of the numbers, and the resolution of raster files
        # we deal with, this shouldn't be a problem.
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

# Function to retrieve coordinate reference from both user and climate data
EC_DataProjection <- function(data, climate.data) {
    if (!is.null(data) & !compareCRS(data, climate.data, verbatim=TRUE)) {
        sp <- SpatialPoints(data)
        if (is.na(crs(sp))) {
            crs(sp) <- '+init=epsg:4326'
        }
        newdata <- as.data.frame(spTransform(sp, crs(climate.data))) # convert to data frame
        names(newdata) <- names(data)
        return(newdata)
    }
    return(data)
}

EC_OutfileName <- function(filename, id_str, ext) {
    return(sprintf("%s_%s.%s", filename, id_str, ext))
}

## Use BIOMOD2 package to format data
# BIOMOD_FormatingData(resp.var, expl.var, resp.xy = NULL, resp.name = NULL, eval.resp.var = NULL,
#   eval.expl.var = NULL, eval.resp.xy = NULL, PA.nb.rep = 0, PA.nb.absences = 1000, PA.strategy = 'random',
#   PA.dist.min = 0, PA.dist.max = NULL, PA.sre.quant = 0.025, PA.table = NULL, na.rm = TRUE)
#
# resp.var a vector, SpatialPointsDataFrame (or SpatialPoints if you work with `only presences' data) containing species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to build the species distribution models.
# expl.var a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to build your models.
# resp.xy optional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to build your models.
# eval.resp.var a vector, SpatialPointsDataFrame your species data (a single species) in binary format (ones for presences, zeros for true absences and NA for indeterminated ) that will be used to evaluate the models with independant data (or past data for instance).
# eval.expl.var a matrix, data.frame, SpatialPointsDataFrame or RasterStack containing your explanatory variables that will be used to evaluate the models with independant data (or past data for instance).
# eval.resp.xy opional 2 columns matrix containing the X and Y coordinates of resp.var (only consider if resp.var is a vector) that will be used to evaluate the modelswith independant data (or past data for instance).
# resp.name response variable name (character). The species name.
# PA.nb.rep number of required Pseudo Absences selection (if needed). 0 by Default.
# PA.nb.absences number of pseudo-absence selected for each repetition (when PA.nb.rep > 0) of the selection (true absences included)
# PA.strategy strategy for selecting the Pseudo Absences (must be `random', `sre', `disk' or `user.defined')
# PA.dist.min minimal distance to presences for `disk' Pseudo Absences selection (in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.dist.max maximal distance to presences for `disk' Pseudo Absences selection(in meters if the explanatory is a not projected raster (+proj=longlat) and in map units (typically also meters) when it is projected or when explanatory variables are stored within table )
# PA.sre.quant quantile used for `sre' Pseudo Absences selection
# PA.table a matrix (or a data.frame) having as many rows than resp.var values. Each column correspund to a Pseudo-absences selection. It contains TRUE or FALSE indicating which values of resp.var will be considered to build models. It must be used with `user.defined' PA.strategy.
# na.rm logical, if TRUE, all points having one or several missing value for environmental data will be removed from analysis


# This uses the biomods function BIOMOD_FormatingData to format user input data.
# It generates pseudo absence points if true absence data are not available or
# adds pseudo absence data to an existing absence dataset.
EC_FormatDataBIOMOD2 <- function(true.absen=NULL,
                                  pseudo.absen.points=0,
                                  pseudo.absen.strategy='random',
                                  pseudo.absen.disk.min=0,
                                  pseudo.absen.disk.max=NULL,
                                  pseudo.absen.sre.quant = 0.025,
                                  climate.data=NULL,
                                  occur=NULL,
                                  species.name=NULL,
                                  save.pseudo.absen=TRUE,
                                  save.env.absen=TRUE,
                                  save.env.occur=TRUE,
                                  generate.background.data=FALSE,
                                  species_algo_str=NULL) {

    # Initialise parameters to default value if not specified
    if (is.null(pseudo.absen.strategy)) {
        pseudo.absen.strategy = 'random'
    }
    if (is.null(pseudo.absen.disk.min)) {
        pseudo.absen.disk.min = 0
    }
    if (is.null(pseudo.absen.sre.quant)) {
        pseudo.absen.sre.quant = 0.025
    }

    if (is.null(true.absen)) { # read true absence point if available.
        absen <- data.frame(lon=numeric(0), lat=numeric(0)) # create an empty data frame for background points
        pseudo.absen.rep = 1 # generate pseudo-absence points
        if (!save.pseudo.absen) {
            pseudo.absen.rep = 0
        }
    }
    else {
        absen <- true.absen
        if (!is.null(climate.data) && nrow(true.absen) > 0) { # ensure true absence dataset is in same projection system as climate.
            absen <- EC_DataProjection(true.absen, climate.data)
        }
        pseudo.absen.rep = 0 # do not generate pseudo absence point when true absence points are available
        pseudo.absen.strategy = 'none'
        pseudo.absen.points = nrow(absen)
        cat("No pseudo absence point is generated.")
    }

    # Generate background data as pseudo absence points
    if (pseudo.absen.strategy != 'none' & generate.background.data) {
        biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))
        myBackgrdData <-
            BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                                 expl.var  = climate.data,
                                 resp.name = species.name,
                                 PA.nb.rep = pseudo.absen.rep,
                                 PA.nb.absences = pseudo.absen.points,
                                 PA.strategy = pseudo.absen.strategy,
                                 PA.dist.min = pseudo.absen.disk.min,
                                 PA.dist.max = pseudo.absen.disk.max,
                                 PA.sre.quant = pseudo.absen.sre.quant)

        # Get background data as absence data
        colnames(myBackgrdData@coord) <- c('lon', 'lat')
        absen <- myBackgrdData@coord[c(which(is.na(myBackgrdData@data.species))), c('lon', 'lat')]

        # Do not generate pseudo absence in next call to BIOMOD_FormatingData
        pseudo.absen.rep = 0
        pseudo.absen.strategy = 'none'
    }

    biomod.data <- rbind(occur[,c("lon", "lat")], absen[,c("lon", "lat")])
    biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))

    myBiomodData <-
        BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                             expl.var  = climate.data,
                             resp.xy   = biomod.data,
                             resp.name = species.name,
                             PA.nb.rep = pseudo.absen.rep,
                             PA.nb.absences = pseudo.absen.points,
                             PA.strategy = pseudo.absen.strategy,
                             PA.dist.min = pseudo.absen.disk.min,
                             PA.dist.max = pseudo.absen.disk.max,
                             PA.sre.quant = pseudo.absen.sre.quant)

    # Save the pseudo absence points generated to file
    pa_filename <- EC_OutfileName(filename="pseudo_absences", id_str=species_algo_str, ext="csv")
    absenv_filename <- EC_OutfileName(filename="absence_environmental", id_str=species_algo_str, ext="csv")
    occenv_filename <- EC_OutfileName(filename="occurrence_environmental", id_str=species_algo_str, ext="csv")
    if (pseudo.absen.rep != 0) {
        pseudoAbsen <- myBiomodData@coord[c(which(is.na(myBiomodData@data.species))), c('lon', 'lat')]
        if (save.pseudo.absen & nrow(pseudoAbsen) > 0) {
            EC_WriteCSV(pseudoAbsen, pa_filename, rownames = FALSE)
        }

        # save the pseudo absence points with environmental variables
        if (save.env.absen) {
            EC_MergeSave(climate.data, pseudoAbsen, species.name, absenv_filename)
        }
    }
    else if (nrow(absen) > 0) {
        # save true-absence/background data generated
        if (!is.null(true.absen)) {
            pa_filename = EC_OutfileName(filename="absence", id_str=species_algo_str, ext="csv") # rename true-absence file
        }
        EC_WriteCSV(absen, pa_filename, rownames = FALSE)

        if (save.env.absen) {
            EC_MergeSave(climate.data, absen, species.name, absenv_filename) # save the true absence points/background points with environmental variables
        }
    }

    if (save.env.occur) {
        EC_MergeSave(climate.data, occur, species.name, occenv_filename) # save the occurrence datasets with environmental variables
    }

    return(myBiomodData)
}

# merge and save  all the data in a csv file
EC_MergeSave <- function(env, csvdata, spname, ofname) {
  data <- cbind(csvdata, species=spname, extract(env, csvdata[c('lon','lat')]))

  EC_WriteCSV(data, ofname, rownames=FALSE)
}

# warning was doing odd things. I just want to print the deng thing.
EC_LogWarning <-function(str, prefix="EcoCommons Warning: ") {
    print(paste(prefix, str, sep=""))
}

EC_RasterLoad <- function(filename) {
    # load raster and assign crs if missing
    r = raster(filename)
    if (is.na(crs(r))) {
        crs(r) = CRS("+init=epsg:4326")
    }
    return(r)
}


# rasters: a vector of rasters, all rasters should have same resolution
# common.crs: crs to use to calculate intersection
EC_RasterExtent <- function(rasters, common.crs) {
    # bring all rasters into common crs (same projection)
    extent.list = lapply(rasters, function(r) { extent(projectExtent(r, common.crs)) })
    # intersect all extents
    common.extent = Reduce(rgeos::intersect, extent.list)
    # compare all against common extents to find out if all extents are the same (used to print warning)
    equal.extents = all(sapply(extent.list, function (x) common.extent == x))

    return (list(equal.extents=equal.extents, common.extent=common.extent))
}

EC_RasterExtentToSTR <- function(ext) {
  return(sprintf("xmin=%f xmax=%f ymin=%f ymax=%f", ext@xmin, ext@xmax, ext@ymin, ext@ymax));
}

# resamplingflag: a flag to determine which resampling approach to take
# selected_layers: a list of indexes to the raster layers to be considered when determine the resolution to be used.
EC_RasterResolution <- function(rasters, resamplingflag, selected_layers) {
    # Return the resolution of the raster given by the index
    get.resolution <- function(i, rasters) {
        return(res(rasters[[i]]))
    }

    # Get resolutions of the raster layers
    if (is.null(selected_layers)) {
        resolutions = lapply(rasters, res)
    }
    else {
        # get the resolutions of the of the selected raster layers only 
        resolutions = lapply(selected_layers, get.resolution, rasters)
    }
    
    if (resamplingflag == "highest") {
        common.res = Reduce(pmin, resolutions)
    } else if (resamplingflag == "lowest") {
        common.res = Reduce(pmax, resolutions)
    }

    # Get resolutions of all input raster layers
    resolutions = lapply(rasters, res)
    is.same.res = all(sapply(resolutions, function(x) all(common.res == x)))
    return (list(common.res=common.res, is.same.res=is.same.res))
}

# generate reference raster with common resolution, crs and extent
EC_RasterRef <- function(rasters, resamplingflag, selected_layers) {
  empty.rasters <- lapply(rasters, function(x) { projectExtent(x, crs(x)) })  # create list of empty rasters to speed up alignment
  common.crs <- crs(empty.rasters[[1]]) # choose a common.crs if all crs in rasters are the same use that one
  
  if (! do.call(compareRaster, c(empty.rasters, extent=FALSE, rowcol=FALSE, prj=TRUE, res=FALSE, orig=FALSE, rotation=FALSE, stopiffalse=FALSE))) {
    common.crs = CRS("+init=epsg:4326") # common prjected adopted in EcoCommons
    EC_LogWarning(sprintf("Auto projecting to common CRS: %s", common.crs))
    empty.rasters = lapply(empty.rasters, function(x) { projectExtent(x, common.crs) })
  }
  
  ce <- EC_RasterExtent(empty.rasters, common.crs)
  if (! ce$equal.extents) {
    EC_LogWarning(sprintf("Auto cropping to common extent %s", EC_RasterExtentToSTR(ce$common.extent)))
  }
  
  cr <- EC_RasterResolution(empty.rasters, resamplingflag, selected_layers)
  if (! cr$is.same.res) {
    EC_LogWarning(sprintf("Auto resampling to %s resolution [%f %f]", resamplingflag, cr$common.res[[1]], cr$common.res[[2]]))
  }
  
  empty.rasters <- lapply(
    empty.rasters,
    function(x) {
      extent(x) = ce$common.extent
      res(x) = cr$common.res
      return(x)
    })
  
  return(empty.rasters[[1]])  # from now on all empty.rasters should be in the same projection
}

EC_RasterWarp <- function(raster.filenames, raster.types, reference, overwrite=TRUE) {
  
  rasters <- mapply(
    function(filename, filetype) {
      gdinfo <- rgdal::GDALinfo(filename)
      mdata <- attr(gdinfo, 'df')
      dtype <- as.character(mdata[['GDType']])
      hasNoDataValues <- mdata[['hasNoDataValue']]
      
      r <- EC_RasterLoad(filename) # warp, crop and rescale raster file if necessary
      dir <- dirname(filename)
      temp_raster <- file.path(dir, paste0(basename(tempfile()), '.tif'))
      te <- extent(reference)
      
      # This is to fix issue with NA value being treated as value 0 if nodatavalue is not set.
      if (hasNoDataValues) {
        gdalUtils::gdalwarp(filename, temp_raster,
                            s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                            te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                            ts=c(ncol(reference), nrow(reference)),
                            # tr=c(...), ... either this or ts
                            r="near",
                            of="GTiff",
                            dstnodata=mdata[['NoDataValue']],
                            co=c("TILED=YES", "COMPRESS=LZW")
        )
      }
      else {
        gdalUtils::gdalwarp(filename, temp_raster,
                            s_srs=CRSargs(crs(r)), t_srs=CRSargs(crs(reference)),
                            te=c(te@xmin, te@ymin, te@xmax, te@ymax),
                            ts=c(ncol(reference), nrow(reference)),
                            # tr=c(...), ... either this or ts
                            r="near",
                            of="GTiff",
                            co=c("TILED=YES", "COMPRESS=LZW")
        )
      }
      
      rasterfilename <- temp_raster
      if (overwrite) {
        file.rename(temp_raster, filename)
        rasterfilename = filename 
      }
      
      r <- raster(rasterfilename)
      if (filetype == "categorical") { # convert to factor if categorical
        r = as.factor(r)
      }
      return(r)
    },
    raster.filenames, raster.types)
  return(rasters)
}

# Load rasters, determine common projection and convert categorical
# rasters to factors
# components:
# raster.filenames is a vector of filenames that will be loaded as rasters
# resamplingflag is a flag to determine which resampling approach to take

EC_RasterResampled  <- function(raster.filenames, raster.types, resamplingflag, selected_layers=NULL, overwrite=TRUE) {
  rasters <- lapply(raster.filenames, EC_RasterLoad)
  
  reference <- EC_RasterRef(rasters, resamplingflag, selected_layers) # determine common raster shape
  
  rasters <- EC_RasterWarp(raster.filenames, raster.types, reference, overwrite) # adjust rasters spatially and convert categorical rasters to factors
  
  return(rasters)
}

# Stack raster input files from user environment
# selected_layers is a list of layers to be considered when determine the resolution of the raster. If none, consider all layers.

EC_EnviroStack <- function(filenames, types, layernames, resamplingflag, selected_layers=NULL) {
  rasters <- EC_RasterResampled (filenames, types, resamplingflag, selected_layers) # adjust rasters to same projection, resolution and extent
  
  rasterstack <- stack(rasters)
  
  names(rasterstack) <- unlist(layernames)  # assign predefined variable names
  return(rasterstack)
}

# Remove raster object and its associated raster files (i.e. grd and gri) if any
EC_RevRasterObject <- function(rasterObject) {
    raster_filenames <- raster_to_filenames(rasterObject, unique = TRUE)
    for (fname in raster_filenames) {
        if (extension(fname)  == '.grd') {
            file.remove(fname, extension(fname, '.gri'))
        }
    }
    rm(rasterObject)
}

EC_SpDataProjection <- function(data, climate.data) {
    sp <- SpatialPoints(data)
    if (is.na(crs(sp))) {
        crs(sp) <- '+init=epsg:4326'
    }

    if (!compareCRS(sp, climate.data, verbatim=TRUE)) {
        sp <- spTransform(sp, crs(climate.data)) # project sp data to the same crs as climate data
    }
    return(sp)
}


# Function to crop the raster to the extent of the constraint region and 
# mask it 

EC_Mask <- function(raster, parsed_geojson) {
  
  cropped_raster <- crop(raster, extent(parsed_geojson))  # crop to the extent of the constraint region before masking
  
  envraster_filename <- paste(EC.env$workdir, basename(tempfile(fileext = ".grd")), sep="/")
  masked_raster <- mask(cropped_raster, parsed_geojson, filename = envraster_filename)
  
  if (is.factor(masked_raster)) {
    masked_raster = as.factor(masked_raster)
  }
  
  EC_RevRasterObject(stack(cropped_raster))  # remove cropped raster and associated raster files (i.e. grd and gri)
  
  return(masked_raster)
}


# Function to geographically constrained the species distribution modelling
# The constraint is the intersection between the occurrence's convex-hull 
# polygon and the constraint. Otherwise, actual constraint is the convex-hull 
# polygon.
# A buffer of width 1-resolution cell is applied to avoid missing environment
# values along the boundary of the polygon.

EC_SDMGeoConstrained <- function(rasterstack, occur, absen, raw_geojson, generateCHull) {
  
  if (is.null(raw_geojson) & !generateCHull) {
    return(list("raster" = rasterstack, "occur" = occur, "absen" = absen))
  }
  
  geojson_crs <- CRS("+init=epsg:3857") # create a geojson for convex-hull polygon if no geojson
  if (is.null(raw_geojson)) {
    parsed_geojson <- SpatialPolygons(list(Polygons(list(Polygon(rbind(c(1,1)))), ID=1)), proj4string=crs(rasterstack))
  } else {
    parsed_geojson <- readOGR(dsn = raw_geojson, layer = "OGRGeoJSON", verbose = FALSE)
    geojson_crs <- crs(parsed_geojson)
  }
  
  if (!compareCRS(rasterstack, parsed_geojson, verbatim=TRUE)) {
    parsed_geojson <- spTransform(parsed_geojson, crs(rasterstack))  # if CRS is different, reproject geojson to rasterstack
  }
  
  if (!is.null(occur)) {
    
    occurSP <- SpatialPoints(occur)  # constrain the occurrence points
    if (is.na(crs(occurSP))) {
      crs(occurSP) <- '+init=epsg:4326'
    }
    
    if (!compareCRS(occurSP, parsed_geojson, verbatim=TRUE)) {
      occurSP <- spTransform(occurSP, crs(parsed_geojson))
    }
    
    region_offset = 0
    if (generateCHull) {
      # get the offset 
      constraint_json <- rjson::fromJSON(raw_geojson)
      region_offset <- constraint_json$properties$region_offset
      if (is.null(region_offset) || is.na(region_offset) || region_offset == '') {
        region_offset = 0
      }
      else {
        region_offset <- as.double(region_offset)
        region_offset <- ifelse(!is.na(region_offset) && is.numeric(region_offset), region_offset/111.0, 0) # convert from km to degree
      }
      
      chcoords <- occurSP@coords[chull(occurSP@coords[,1:2]),]
      chullPolygon <- SpatialPolygons(list(Polygons(list(Polygon(chcoords[,1:2], hole=FALSE)), ID=1)), proj4string=crs(parsed_geojson))
      if (!is.null(raw_geojson)) {
        parsed_geojson <- intersect(parsed_geojson, chullPolygon)
      }
      else {
        parsed_geojson <- chullPolygon
      }
    }
    
    parsed_geojson <- gBuffer(parsed_geojson, byid=TRUE, width=max(region_offset, max(res(rasterstack@layers[[1]])))) #create a buffer of width 1-resolution cell
    
    if (generateCHull) {
      filename <- file.path(EC.env$outputdir, 'modelling_region.json')
      transformed_geojson <- spTransform(parsed_geojson, geojson_crs) 
      writeOGR(transformed_geojson, filename, 'OGRGeoJSON', driver='GeoJSON')  # save the convex-hull generated as geojson.
      transformed_geojson = NULL
      
      gjson <- rjson::fromJSON(file=filename)
      origjson <- rjson::fromJSON(raw_geojson)
      if (! 'crs' %in% names(gjson)) {
        gjson = append(gjson, list(crs=origjson$crs))   
      }
      if (! 'properties' %in% names(gjson)) {
        gjson = append(gjson, list(properties=origjson$properties))
      }
      write(rjson::toJSON(gjson), filename)
      
      origjson = NULL
      gjson = NULL
    }    
    
    occurSPconstrained <- occurSP[!is.na(over(occurSP, parsed_geojson))]
    occurconstrained <- as.data.frame(occurSPconstrained)
    
    absenconstrained = NULL
    
    if (!is.null(absen) && nrow(absen) > 0) {
      absenSP <- EC_SpDataProjection(absen, rasterstack)
      absenSPconstrained <- absenSP[!is.na(over(absenSP, parsed_geojson))]
      
      absenSPconstrained <- spTransform(absenSPconstrained, CRS('+init=epsg:4326')) # project it back to epsg:4326 for saving as a csv file
      absenconstrained <- as.data.frame(absenSPconstrained)
      #names(absenconstrained) <- c("lon", "lat")
      #write.csv(absen, file=absenFilename, row.names=FALSE)
    }
  }
  else {
    occurconstrained = NULL
    absenconstrained = NULL
  }
  
  geoconstrained <- stack(lapply(as.list(rasterstack), EC_Mask, parsed_geojson))  # mask the rasterstack
  
  mylist <- list("raster" = geoconstrained, "occur" = occurconstrained, "absen" = absenconstrained)
  return(mylist)
}

#------------------------------------------------------------------------------------------------------------
EC_CreateBackgroundData <- function(rastermask, bgsize) {
    background <- dismo::randomPoints(rastermask, n=bgsize, extf=1.0, excludep=FALSE)
    colnames(background) <- c('lon', 'lat')
    return(as.data.frame(background))
}

# If chosen to use a bias file, apply it in the occurences to the background
# points 

EC_CreateBiasedBackgroundData <- function(bg, biasfile, bgsize) {
  
  bias <- raster(biasfile)
  
  bg_sp_class<-SpatialPoints(coords=cbind(bg$lon,bg$lat), proj4string=CRS(as.character(NA)))   # convert background points to spatial points formal class
  
  weight <- raster::extract(bias, bg_sp_class)  
  bg_sp_class <- bg_sp_class[!is.na(weight), ]
  weight <- weight[!is.na(weight)] # random sample of the candidate background points
  bg_biased <- bg_sp_class[sample(length(bg_sp_class), size=bgsize, replace = TRUE, prob=weight), ]
  colnames(bg_biased@coords) <- c("lon", "lat")
  bg_biased <- bg_biased@coords
  return(bg_biased)
}

# Plot a response curve to show results

EC_PlotResponseCurve <- function(model, predvars, species, oformat) {  # overall response curves
  png(filename <- file.path(EC.env$outputdir, paste0("response_curves_", species, "_maxent_", oformat, ".png")),
      width=800, height=800)
  response(model, at=mean, range='p', expand=5) 
  dev.off()
  
  for (varname in predvars) {  # save individual response curve
    png(filename=file.path(EC.env$outputdir, paste0("response_curve_", species, "_maxent_", oformat, "_", varname, ".png")),
        width=800, height=800)
    response(model, var=varname, at=mean, 
             fun=function(x, y, ...) predict(x, y, args=c(outputformat=oformat)), 
             range='p', expand=5)
    dev.off()    
  }
}

#------------------------------------------------------------------------------------------------------------

# Load and export plot function for continuous predictor variables for 
# the occurrence data

EC_PlotContinOCC <- function(dataToExplore, listVars, doCorrelation=TRUE, 
                             outerTitle="", fnamePrefix="occurrence_predictors_", outputdir="output")  {
  
  nVars = length(listVars)
  
  if (nVars == 0) {
    return()
  }
  
  nCols = ceiling(sqrt(nVars))
  nRows = ceiling(nVars/nCols)
  plotWidth = 800; plotHeight = 800
  
  # Boxplots
  png(filename=file.path(outputdir, paste0(fnamePrefix, 'boxplot.png')),
      width=plotWidth, height=plotHeight)
  par(mfrow=c(nRows,nCols), oma=c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    boxplot(dataToExplore[,c(nameVar)], main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Histograms
  png(filename=file.path(outputdir, paste0(fnamePrefix, 'histogram.png')),
      width=plotWidth, height=plotHeight)
  par(mfrow=c(nRows,nCols), oma=c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    hist(dataToExplore[,c(nameVar)], breaks=10, xlab="", main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Density plots
  png(filename=file.path(outputdir, paste0(fnamePrefix, 'density.png')),
      width=plotWidth, height=plotHeight)
  par(mfrow=c(nRows,nCols), oma=c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    plot(density(dataToExplore[,c(nameVar)], na.rm=TRUE), main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Pearson Correlation
  if (doCorrelation && nCols > 1) {
    
    png(filename=file.path(outputdir, paste0(fnamePrefix, 'correlations.png')),
        width=plotWidth, height=plotHeight)
    par(oma=c(0,0,1.5,0), ps=15)
    pairs(dataToExplore, 
          lower.panel = function(x,y){
            usrSaved = par("usr"); on.exit(par(usrSaved))
            par(usr = c(0, 1, 0, 1), ps=15)
            text(0.5, 0.5, cex=2,
                 round(cor(x, y, use="na.or.complete", method="pearson"), digits=2))})
    title(outerTitle, outer=TRUE, line=0)
    dev.off()
  }
}  

# Load and export plot function for categorical predictor variables for 
# the occurrence data

EC_PlotCategOCC <- function(dataToExplore, listVars, outerTitle="", 
                            fnamePrefix="occurrence_predictors_", outputdir="output") {
  
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
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    barplot(table(dataToExplore[,c(nameVar)]), 
            main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
}


# Load and export plot function for continuous predictor variables for 
# the background data

EC_PlotContinBG <- function(dataToExplore, listVars, doCorrelation=TRUE, outerTitle="", fnamePrefix="background_predictors_", outputdir="output")  {
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
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    boxplot(dataToExplore[,c(nameVar)], main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Histograms
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'histogram.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    hist(dataToExplore[,c(nameVar)], breaks=10, xlab="", main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
  
  # Density plots
  png(filename = file.path(outputdir, paste0(fnamePrefix, 'density.png')),
      width = plotWidth, height = plotHeight)
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    plot(density(dataToExplore[,c(nameVar)], na.rm=TRUE), main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
}

# Load plot function for categorical predictor variables for the background data
EC_PlotCategBG <- function(dataToExplore, listVars, outerTitle="",
                           fnamePrefix="background_predictors_", outputdir="output") {

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
  par(mfrow = c(nRows,nCols), oma = c(0,0,1.5,0), ps=15)
  for (nameVar in listVars) {
    barplot(table(dataToExplore[,c(nameVar)]), 
            main=nameVar)
  }
  title(outerTitle, outer=TRUE, line=0)
  dev.off()
}

#------------------------------------------------------------------------------------------------------------
# Summary statistics to evaluate model performance
# Return equal sensitivity specificity threshold

EC_SummaryStats <- function(occur, absen, model, predictor, species) {
    writeLines('Compute evaluation statistics ...')
    # Get evaluation statistics for different threshold methods
    evaluation<-dismo::evaluate(occur, absen, model, predictor)

    png(file=file.path(EC.env$outputdir, paste0("density_", species, "_maxent.png")),
        width=800, height=800)
    density(evaluation)
    dev.off()

    # Retrieve threshold values for each method
    threshold_all <- threshold(evaluation)
    threshold_all <- subset(threshold_all[,1:5]) # leave out sensitivity threshold
    threshold_kappa <- threshold(evaluation, stat='kappa')
    threshold_spec_sens <- threshold(evaluation, stat='spec_sens')
    threshold_no_omission <- threshold(evaluation, stat='no_omission')
    threshold_prevalence <- threshold(evaluation, stat='prevalence')
    threshold_equal_sens_spec <- threshold(evaluation, stat='equal_sens_spec')

    # Run the evaluate function for each threshold method
    evaluation_kappa <- dismo::evaluate(occur, absen, model, predictor, threshold_kappa)
    evaluation_spec_sens <- dismo::evaluate(occur, absen, model, predictor, threshold_spec_sens)
    evaluation_no_omission <- dismo::evaluate(occur, absen, model, predictor, threshold_no_omission)
    evaluation_prevalence <- dismo::evaluate(occur, absen, model, predictor, threshold_prevalence)
    evaluation_equal_sens_spec <- dismo::evaluate(occur, absen, model, predictor, threshold_equal_sens_spec)

    # Create a dataframe of evaluation statistics using the threshold where Kappa is highest
    kappa_prevalence = t(as.data.frame(as.double(slot(evaluation_kappa,"prevalence"))))
    kappa_OR = t(as.data.frame(as.double(slot(evaluation_kappa,"OR"))))
    kappa_MCR = t(as.data.frame(as.double(slot(evaluation_kappa,"MCR"))))
    kappa_NPP = t(as.data.frame(as.double(slot(evaluation_kappa,"NPP"))))
    kappa_PPP = t(as.data.frame(as.double(slot(evaluation_kappa,"PPP"))))
    kappa_FNR = t(as.data.frame(as.double(slot(evaluation_kappa,"FNR"))))
    kappa_FPR = t(as.data.frame(as.double(slot(evaluation_kappa,"FPR"))))
    kappa_TNR = t(as.data.frame(as.double(slot(evaluation_kappa,"TNR"))))
    kappa_TPR = t(as.data.frame(as.double(slot(evaluation_kappa,"TPR"))))
    kappa_CCR = t(as.data.frame(as.double(slot(evaluation_kappa,"CCR"))))
    kappa_ODP = t(as.data.frame(as.double(slot(evaluation_kappa,"ODP"))))

    kappa_bind <- rbind(kappa_CCR, kappa_FNR, kappa_FPR, kappa_MCR, kappa_NPP, kappa_ODP, kappa_OR, kappa_PPP, kappa_prevalence, kappa_TNR, kappa_TPR)


    # Create a dataframe of evaluation statistics using the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
    spec_sens_prevalence = t(as.data.frame(as.double(slot(evaluation_spec_sens,"prevalence"))))
    spec_sens_OR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"OR"))))
    spec_sens_MCR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"MCR"))))
    spec_sens_NPP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"NPP"))))
    spec_sens_PPP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"PPP"))))
    spec_sens_FNR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"FNR"))))
    spec_sens_FPR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"FPR"))))
    spec_sens_TNR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"TNR"))))
    spec_sens_TPR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"TPR"))))
    spec_sens_CCR = t(as.data.frame(as.double(slot(evaluation_spec_sens,"CCR"))))
    spec_sens_ODP = t(as.data.frame(as.double(slot(evaluation_spec_sens,"ODP"))))

    spec_sens_bind <- rbind(spec_sens_CCR, spec_sens_FNR, spec_sens_FPR, spec_sens_MCR, spec_sens_NPP, spec_sens_ODP, spec_sens_OR, spec_sens_PPP, spec_sens_prevalence, spec_sens_TNR, spec_sens_TPR)

    # Create a dataframe of evaluation statistics using the highest threshold at which there is no omission
    no_omission_prevalence = t(as.data.frame(as.double(slot(evaluation_no_omission,"prevalence"))))
    no_omission_OR = t(as.data.frame(as.double(slot(evaluation_no_omission,"OR"))))
    no_omission_MCR = t(as.data.frame(as.double(slot(evaluation_no_omission,"MCR"))))
    no_omission_NPP = t(as.data.frame(as.double(slot(evaluation_no_omission,"NPP"))))
    no_omission_PPP = t(as.data.frame(as.double(slot(evaluation_no_omission,"PPP"))))
    no_omission_FNR = t(as.data.frame(as.double(slot(evaluation_no_omission,"FNR"))))
    no_omission_FPR = t(as.data.frame(as.double(slot(evaluation_no_omission,"FPR"))))
    no_omission_TNR = t(as.data.frame(as.double(slot(evaluation_no_omission,"TNR"))))
    no_omission_TPR = t(as.data.frame(as.double(slot(evaluation_no_omission,"TPR"))))
    no_omission_CCR = t(as.data.frame(as.double(slot(evaluation_no_omission,"CCR"))))
    no_omission_ODP = t(as.data.frame(as.double(slot(evaluation_no_omission,"ODP"))))

    no_omission_bind <- rbind(no_omission_CCR, no_omission_FNR, no_omission_FPR, no_omission_MCR, no_omission_NPP, no_omission_ODP, no_omission_OR, no_omission_PPP, no_omission_prevalence, no_omission_TNR, no_omission_TPR)

    # Create a dataframe of evaluation statistics using the threshold at which the modeled prevalence is closest to the observed prevalence
    prevalence_prevalence = t(as.data.frame(as.double(slot(evaluation_prevalence,"prevalence"))))
    prevalence_OR = t(as.data.frame(as.double(slot(evaluation_prevalence,"OR"))))
    prevalence_MCR = t(as.data.frame(as.double(slot(evaluation_prevalence,"MCR"))))
    prevalence_NPP = t(as.data.frame(as.double(slot(evaluation_prevalence,"NPP"))))
    prevalence_PPP = t(as.data.frame(as.double(slot(evaluation_prevalence,"PPP"))))
    prevalence_FNR = t(as.data.frame(as.double(slot(evaluation_prevalence,"FNR"))))
    prevalence_FPR = t(as.data.frame(as.double(slot(evaluation_prevalence,"FPR"))))
    prevalence_TNR = t(as.data.frame(as.double(slot(evaluation_prevalence,"TNR"))))
    prevalence_TPR = t(as.data.frame(as.double(slot(evaluation_prevalence,"TPR"))))
    prevalence_CCR = t(as.data.frame(as.double(slot(evaluation_prevalence,"CCR"))))
    prevalence_ODP = t(as.data.frame(as.double(slot(evaluation_prevalence,"ODP"))))

    prevalence_bind <- rbind(prevalence_CCR, prevalence_FNR, prevalence_FPR, prevalence_MCR, prevalence_NPP, prevalence_ODP, prevalence_OR, prevalence_PPP, prevalence_prevalence, prevalence_TNR, prevalence_TPR)

    # Create a dataframe of evaluation statistics using the threshold where the sensitivity and specificity are equal
    equal_sens_spec_prevalence = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"prevalence"))))
    equal_sens_spec_OR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"OR"))))
    equal_sens_spec_MCR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"MCR"))))
    equal_sens_spec_NPP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"NPP"))))
    equal_sens_spec_PPP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"PPP"))))
    equal_sens_spec_FNR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"FNR"))))
    equal_sens_spec_FPR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"FPR"))))
    equal_sens_spec_TNR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"TNR"))))
    equal_sens_spec_TPR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"TPR"))))
    equal_sens_spec_CCR = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"CCR"))))
    equal_sens_spec_ODP = t(as.data.frame(as.double(slot(evaluation_equal_sens_spec,"ODP"))))

    equal_sens_spec_bind <- rbind(equal_sens_spec_CCR, equal_sens_spec_FNR, equal_sens_spec_FPR, equal_sens_spec_MCR, equal_sens_spec_NPP, equal_sens_spec_ODP, equal_sens_spec_OR, equal_sens_spec_PPP, equal_sens_spec_prevalence, equal_sens_spec_TNR, equal_sens_spec_TPR)

    # Generate evaluation statistics table
    all_thresholds_matrix <- cbind(kappa_bind, spec_sens_bind, no_omission_bind, prevalence_bind, equal_sens_spec_bind)

    # Add in the original threshold values
    colnames(all_thresholds_matrix) <- colnames(threshold_all)
    threshold_matrix <- rbind(threshold_all, all_thresholds_matrix) 
    
    # Change column and row names
    colnames(threshold_matrix) <- c("Max Kappa", "Max sum sensitivity and specificity", "No Omission", "Prevalence", "Equal sensitivity specificity") 
    rownames(threshold_matrix) <- c("Threshold value", "Prevalence","Odds ratio","Misclassification rate","Negative predictive power","Positive predictive power","False negative rate","False positive rate","True negative rate","True positive rate","Correct classification rate","Overall diagnostic power")

    # Write evaluation statistics as csv
    write.csv(threshold_matrix, file=file.path(EC.env$outputdir, paste0('evaluation_', species, '_maxent.csv')))

    # return equal sensitivity specificity threshold
    return(threshold_equal_sens_spec)
}


#------------------------------------------------------------------------------------------------------------
## Computes the Area/Extent of Occuppancy (AOO/EOO) for species of interest

EC_AreaExtentOccupancy <- function(occur, species) {
  writeLines('Compute Area Extent Occuopancy ...')
  data(land) #land is built in to the package, it is possible to add own shapefile if needed
  
  cnames = names(occur)
  occur$Species = species
  if (!('year' %in% cnames)) {
    occur = occur[, c('lat','lon','species')]
    
    # Get AOO value
    aoo <- AOO.computing(occur, Cell_size_AOO=2, nbe.rep.rast.AOO=0,
                         parallel=FALSE, NbeCores =2, show_progress=FALSE, export_shp=FALSE)
    
    # Generate table with AOO and EOO values
    aoo_eoo <- IUCN.eval(occur, file_name=file.path(EC.env$outputdir, paste0("aoo_eoo_", species, ".csv")))
    aoo_eoo <- aoo_eoo[c('taxa', 'EOO', 'AOO', 'Nbe_unique_occ.')]
    colnames(aoo_eoo) <- c('taxa', 'EOO', 'AOO', 'Number of unique occurrences')
    write.csv(aoo_eoo, file=file.path(EC.env$outputdir, paste0("aoo_eoo_", species, "_maxent.csv")))
    
    # Run IUCN evaluation to generate map
    IUCN.eval(DATA = occur, country_map = land, Cell_size_AOO = 2, Cell_size_locations = 10,
              Resol_sub_pop = 5, method_locations = "fixed_grid", Rel_cell_size = 0.05,
              DrawMap = TRUE, add.legend = TRUE,
              file_name = NULL, export_shp = FALSE, write_shp = FALSE,
              write_results=TRUE, protec.areas = NULL, map_pdf = FALSE, draw.poly.EOO=TRUE,
              exclude.area = FALSE, method_protected_area = "no_more_than_one",
              buff_width = 0.1, SubPop=FALSE, alpha=1, buff.alpha=0.1,
              method.range="convex.hull", nbe.rep.rast.AOO=0,
              showWarnings=TRUE, write_file_option="csv",
              parallel=FALSE, NbeCores=2)
  }
  else {
    occur = occur[, c('lat', 'lon', 'species', 'year')] # AOO and EOO evaluation per year, if applicable
    
    # Subsetting data.frames into every 5 years until the last ten years which is per year
    yrmin = min(occur['year'])
    yrmax = max(occur['year'])
    if ((yrmax - yrmin) < 10) {
      ylist = seq(yrmin, yrmax, 1)
    } else {
      ylist = seq(yrmax-9, yrmax, 1)
    }
    sub_list = lapply(ylist, function(x) occur[occur$year == x,]) 
    ynames = lapply(ylist, function(x) paste0(x))
    names(sub_list) <- ynames
    
    # 5-year group
    if ((yrmax - yrmin) >= 10) {  
      remaining = (yrmax-yrmin) %% 5
      yrmin = yrmin - (4 - remaining)
      ylist = seq(yrmin, yrmax-10, 5)
      sub_list2 = lapply(ylist, function(x) occur[occur$year>=x & occur$year <= (x+4),])
      ynames = lapply(ylist, function(x) paste0(x, '-', x+4))
      names(sub_list2) <- ynames
      sub_list = c(sub_list2, sub_list)
    }
    
    # Select only those years with data, and with 3 or more occurrence records 
    # (i.e. minimum required to get EOO result)
    subset <- sub_list[lapply(sub_list, function(x) dim(x)[1]) > 2]
    
    # Get the AOO, EOO, and IUCN evaluation for each year range and all results in a list
    aoo_eoo_result <- lapply(subset, function(x){
      IUCN.eval(DATA = x, country_map = land, Cell_size_AOO = 2, Cell_size_locations = 10,
                Resol_sub_pop = 5, method_locations = "fixed_grid", Rel_cell_size = 0.05,
                DrawMap = TRUE, add.legend = TRUE,
                file_name = NULL, export_shp = FALSE, write_shp = FALSE,
                write_results=TRUE, protec.areas = NULL, map_pdf = FALSE, draw.poly.EOO=TRUE,
                exclude.area = FALSE, method_protected_area = "no_more_than_one",
                buff_width = 0.1, SubPop=FALSE, alpha=1, buff.alpha=0.1,
                method.range="convex.hull", nbe.rep.rast.AOO=0,
                showWarnings=TRUE, write_file_option="csv", 
                parallel=FALSE, NbeCores=2) 
    })
    
    aoo_eoo_year <- bind_rows(aoo_eoo_result)  # combine the rows of the list into a data frame
    aoo_eoo_year['year'] = names(aoo_eoo_result)
    
    aoo_eoo_year = aoo_eoo_year[c('year', 'EOO', 'AOO', 'Nbe_unique_occ.')]
    names(aoo_eoo_year) <- c('year', 'EOO', 'AOO', 'Number of unique occurrences')
    ofname = file.path(EC.env$outputdir, paste0('aoo_eoo_year_', species, '_maxtent.csv'))
    write.csv(aoo_eoo_year, ofname, row.names = TRUE)  # the csv should show the year (subset) in the first column and then add the 4 columns ID, EOO, AOO, Nbe_unique_occ
    
    # Plot time series to show trend in AOO and EOO over time
    aoo_time_series = t(data.matrix(aoo_eoo_year['AOO']))
    colnames(aoo_time_series) <- unlist(aoo_eoo_year['year'])
    png(filename=file.path(EC.env$outputdir, paste0('aoo_year_', species, '_maxent.png')),
        width=800, height=800)
    barplot(aoo_time_series, main="Area of Occupancy over time",
            xlab="Year range ", ylab="Area of Occupancy (km2)", pch=19, las=2)
    dev.off()
    
    eoo_time_series = t(data.matrix(aoo_eoo_year['EOO']))
    colnames(eoo_time_series) <- unlist(aoo_eoo_year['year'])
    png(filename=file.path(EC.env$outputdir, paste0('eoo_year_', species, '_maxent.png')),
        width=800, height=800)
    barplot(eoo_time_series, main="Extent of Occupancy over time",
            xlab="Year range ", ylab="Extent of Occupancy (km2)", pch=19, las=2)
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# Generate a map that converts the continuous prediction into a presence/absence
# map based on 'equal sensitivity specificity' threshold. 
# presence-1 (red) >= threshold, absence-0 (blue) < threshold

EC_ThresholdMap <- function(proj_raster, threshold) {
  writeLines('Generating presence/absence map ...')
  present_map <- calc(proj_raster, fun=function(x) { return(as.integer(x >= threshold))})
  
  threshold_map = file.path(EC.env$outputdir, paste0("equal.sens.spec_", occur.species, "_maxent_cloglog.tif"))
  writeRaster(present_map, threshold_map, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  
  threshold_map = file.path(EC.env$outputdir, paste0("equal.sens.spec_", occur.species, "_maxent_cloglog.png"))
  png(filename=threshold_map, width=800, height=800)
  plot(present_map, xlab='longitude', ylab='latitude')
  dev.off()
}

#------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------------------------------------------

# function to plot projection tiff file (with histogram)
EC_PlotProjection <- function(inputfile, main) {
    ## define the breaks of the color key
    my.at <- seq(0,1.0,by=0.1)
    ## the labels will be placed vertically centered
    my.labs.at <- seq(0,1.0,by=0.25)
    ## define the labels
    my.lab <- seq(0,1.0,by=0.25)
    ## define colors
    my.col <- colorRampPalette(c("grey90","yellow4","green4"))(100)

    # Read in tiff input file as rasterstack and plot it
    require('rasterVis')
    levelplot(stack(raster(inputfile)),
              at=my.at,
              margin=T,
              col.regions=my.col,
              main=main,
              colorkey=list(labels=list(
                labels=my.lab,
                at=my.labs.at)))
}

# function to generate a filename for the specified file type and extension.
EC_FilePath <- function(file_type, projection_name, species, 
                        outputdir=EC.env$outputdir, filename_ext=NULL, file_ext='tif') {
  if (is.null(filename_ext)) {
    basename = paste(file_type, projection_name, species, sep="_")
  }  else {
    basename = paste(file_type, projection_name, species, filename_ext, sep="_")
  }
  return(file.path(outputdir, paste(basename, file_ext, sep=".")))
}

EC_ProjectionImage <- function(inputfile, projection.name, species, 
                               species_algo_str, outputdir=EC.env$outputdir, filename_ext=NULL) {
  filename = EC_FilePath("proj", projection.name, species_algo_str, 
                         outputdir, filename_ext, "png")
  png(filename)
  title = paste(species, projection.name, "projections", sep=" ")
  plot(raster(inputfile), main=title, xlab='longitude', ylab='latitude')
  # TODO: to use levelplot to produce histogram instead of plot.
  #EC_PlotProjection(inputfile, title)
  dev.off()
}

# function to compute and save occurrence probability change metrics as geotif file
EC_OccurenceProbMetric <- function(prob_rasters, outfilename) {
    changeproj <- overlay(prob_rasters[[1]], prob_rasters[[2]], fun=function(r1, r2) { return(r1-r2) })
    writeRaster(changeproj, outfilename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
}

# function to compute and save species range change metric as geotif file
EC_SpRangeMetric <- function(prob_rasters, threshold, outfilename) {
    # return 1 for Blank, 3 for Expansion, 0 for Contraction and 2 for No Change
    rangeChange <- overlay(as.integer(prob_rasters[[1]] >= threshold),
                           as.integer(prob_rasters[[2]] >= threshold),
                           fun=function(fp, cp) { return((2 * fp) + 1 - cp)})
    writeRaster(rangeChange, outfilename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)

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

# function to compute and save Centre of Gravity as csv file.
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
    results[1,] = EC_Distance(current_cog['COGy'], current_cog['COGx'], future_cog['COGy'], future_cog['COGx'])
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


# function to save projection output raster
EC_SaveModelProjection <- function(model.obj, projection.name, species, 
                                   species_algo_str, outputdir=EC.env$outputdir, filename_ext=NULL) {
  
  filename = EC_FilePath("proj", projection.name, species_algo_str, outputdir, filename_ext, "tif")
  writeRaster(model.obj, filename, format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
  pngfilename = EC_FilePath("proj", projection.name, species_algo_str, outputdir, filename_ext, "png")
  png(pngfilename)
  title = paste(species, projection.name, "projections", sep=" ")
  plot(model.obj, xlab="latitude", ylab="longtitude", main=title)
  dev.off()
  
  return (filename)
}

# function to save RData in outputdir
EC_Save <- function(robj, name, outputdir=EC.env$outputdir) {
  filename = file.path(outputdir, name)
  save(robj, file=filename)
}

# function to save CSV Data in outputdir
EC_WriteCSV <- function(robj, name, outputdir=EC.env$outputdir, rownames=TRUE) {
    filename = file.path(outputdir, name)
    write.csv(robj, file=filename, row.names=rownames)
}

# function to get model object
EC_GetModelObject <- function(model.file=EC.env$inputmodel) {
  return (get(load(file=model.file)))
}

# convert all .gri/.grd found in folder to gtiff
# TODO: extend to handle other grid file formats, e.g. .asc
EC_GRDtoGTIFF <- function(folder, algorithm, filename_ext=NULL, noDataValue=NULL) {
    grdfiles <- list.files(path=folder,
                           pattern="^.*\\.grd")
    for (grdfile in grdfiles) {
        ext = file_ext(grdfile)
        if (!is.null(ext)) {
            pattern = paste0('\\.', ext, '$')
            grdname <- sub(pattern, '', grdfile)
        }

        grd <- raster(file.path(folder, grdfile))

        if (is.na(proj4string(grd))) {
            crs = CRS("+init=epsg:4326")  # if projection is missing, initialise it to EPSG:4326
            proj4string(grd) <- crs
        }
        
        basename = paste(grdname, algorithm, sep="_")
        if (!is.null(filename_ext)) {
            basename = paste(grdname, algorithm, filename_ext, sep="_")
        }
        filename = file.path(folder, paste(basename, 'tif', sep="."))  # write raster as geotiff

        # To do: This is a temporary fix for nodatavalue is not recognised by mosaic_raster
        # due to a bug in gdal libarry. It shall be removed when using gdal 2.1.3.
        dtype = dataType(grd)
        if (is.null(noDataValue)) {
            writeRaster(grd, filename, datatype=dataType(grd),
                        format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
        }
        else {
            writeRaster(grd, filename, datatype=dataType(grd), NAflag=noDataValue,
                        format="GTiff", options=c("COMPRESS=LZW", "TILED=YES"), overwrite=TRUE)
        }
        # remove grd files
        file.remove(file.path(folder, paste(grdname, c("grd","gri"), sep=".")))
    }
}

############################################################
#
# define helper functions for projections
#
############################################################

# function to check that the environmental layers used to project the
# model are the same as the ones used to create the model object
#    model.obj     ... model to project
#    climatelayers ... climate data to project onto
EC_CheckLayers <- function(model.obj, climatelayers, climate_filenames) {
    message("Checking environmental layers used for projection")
  
    if (inherits(model.obj, "DistModel")) {  # dismo package
        model.layers = colnames(model.obj@presence)
    } else if (inherits(model.obj, "gbm")) {  # brt package
        model.layers = summary(model.obj)$var
    } else if (inherits(model.obj, "BIOMOD.models.out")) {  # biomod package
        model.layers = model.obj@expl.var.names
    }
    pred.layers = names(climatelayers)

    # check if the env layers were in the original model
    if (sum(!(pred.layers %in% model.layers)) > 0 ){
        # To do: Shall remove this sometimes later.
        # The model layer name is to be matched to climate layer name, or its file name.
        # This is to allow for old SDM result to be used.
        if (sum(!(model.layers %in% pred.layers)) > 0){
            filenames = lapply(climate_filenames, function(x) sub("^([^.]*).*", "\\1", basename(x)))
            indexes = match(model.layers, filenames)
            for (i in indexes){
                if (!is.na(i)){
                    pred.layers[i] = model.layers[i]    #Use the corresponding layer name in the model
                }
            }
            names(climatelayers) = pred.layers
        }

        message("Dropping environmental layers not used in the original model creation...")
        # create a new list of env predictors by dropping layers not in the original model
        new.predictors = climatelayers
        for (pl in pred.layers) {
            if (!(pl %in% model.layers)) {
                new.predictors = dropLayer(new.predictors, pl)
            }
        }
        return(new.predictors)
    } else {
        return(climatelayers)
    }
}


EC_FamilyFromString <- function(s)
{
    # get family from a string (character) in a safe way
    # works for all variants of the R family object (e.g. see ?family)
    # i.e.
    # EC_FamilyFromString("binomial")
    # EC_FamilyFromString("binomial(link=logit)")
    # EC_FamilyFromString("binomial(link=\"logit\")")
    # ...
    # EC_FamilyFromString("quasi(link = \"identity\", variance = \"constant\")")

    s=gsub(pattern="\"|| ", replacement="", s) # strip quotes and spaces
    f=gsub(pattern="\\(.*\\)", replacement="", s) # the name of the function

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
               replacement="",
               sub(pattern=f,
                    replacement="",
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

#' Grid Information from Geographic (lat lon) Projections
#'
#' Since spatial grids in geographic projections do not have equal area or
#' perimeters, \code{EC_GridInfo} extracts perimeter & area related information
#' for latitudinal bands with differing longitudinal widths. \cr\cr Outputs
#' lengths are in m using Vincenty's equation (\code{distance})and areas in m2.
#' Surface areas are calculated summing surface areas of spherical polygons as
#' estimated using l'Huiller's formula.
#'
#'
#' @param lats is a vector of latitudes representing the midpoint of grid cells
#' @param cellsize is a single value (assuming square cells) or a two value
#' vector (rectangular cells) representing the height (latitude) and width
#' (longitude) of the cells
#' @param r is a single value representing the radius of the globe in m.
#' Default is for the WGS84 elipsoid
#' @return a data.frame listing: \item{lat}{the latitude representing the
#' midpoint of the cell} \item{top}{length of the top of the cell (m)}
#' \item{bottom}{length of the bottom of the cell (m)} \item{side}{length of
#' the side of the cell (m)} \item{diagnal}{length of the diagnals of the cell
#' (m)} \item{area}{area of the cell (m2)}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @references information on l'Huiller's formula
#' \url{http://williams.best.vwh.net/avform.htm for more info)} code for
#' estimating area of polygon on sphere was modified from
#' \url{http://forum.worldwindcentral.com/showthread.php?t=20724}
#' @examples
#'
#' #show output for latitudes from -87.5 to 87.5 at 5 degree intervals
#' EC_GridInfo(lats=seq(-87.5,87.5,5), 5)
#'
#' @export
# This is a fix for SDM tool EC_GridInfo function due to floating point operation.
EC_GridInfo <- function(lats,cellsize,r=6378137) {
    r2 = r^2 #radius of earth
    ###need checks to ensure lats will not go beyond 90 & -90
    if (length(cellsize)==1) cellsize=rep(cellsize,2) #ensure cellsize is defined for both lat & lon
    out = data.frame(lat=lats) #setup the output dataframe
    toplats = lats+(0.5*cellsize[1]); bottomlats = lats-(0.5*cellsize[1]) #define the top and bottom lats
    check = range(c(toplats,bottomlats),na.rm=TRUE); if (-90.0001>check[1] | 90.0001<check[2]) stop('latitudes must be between -90 & 90 inclusively')
    out$top = EC_Distance(toplats,rep(0,length(lats)),toplats,rep(cellsize[2],length(lats)))$distance
    out$bottom = EC_Distance(bottomlats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
    out$side = EC_Distance(toplats,rep(0,length(lats)),bottomlats,rep(0,length(lats)))$distance
    out$diagnal = EC_Distance(toplats,rep(0,length(lats)),bottomlats,rep(cellsize[2],length(lats)))$distance
    #calculate area of a spherical triangle using spherical excess associated by knowing distances
    #tan(E/4) = sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2))
    #where a, b, c = sides of spherical triangle
    #s = (a + b + c)/2
    #from CRC Standard Mathematical Tables
    #calculate excess based on  l'Huiller's formula (http://williams.best.vwh.net/avform.htm for more info)
    #code modified from (http://forum.worldwindcentral.com/showthread.php?t=20724)
    excess = function(lam1,lam2,beta1,beta2){ #calculate excess... inputs are in radians
        haversine = function(y) { (1-cos(y))/2 }
        cosB1 = cos(beta1); cosB2 = cos(beta2)
        hav1 = haversine(beta2-beta1) + cosB1*cosB2*haversine(lam2-lam1)
        aa = 2 * asin(sqrt(hav1)); bb = 0.5*pi - beta2; cc = 0.5*pi - beta1
        ss = 0.5*(aa+bb+cc)
        tt = tan(ss/2)*tan((ss-aa)/2)*tan((ss-bb)/2)*tan((ss-cc)/2)
        return(abs(4*atan(sqrt(abs(tt)))))
    }
    if (any(bottomlats==-90)) { pos = which(bottomlats==-90); bottomlats[pos] = -bottomlats[pos]; toplats[pos] = -toplats[pos]} #ensure no -90 bottom lats
    out$area = excess(lam1=0,lam2=cellsize[2]*pi/180,toplats*pi/180,toplats*pi/180)
    out$area = abs(out$area-excess(lam1=0,lam2=cellsize[2]*pi/180,bottomlats*pi/180,bottomlats*pi/180))*r2
    return(out)
}

#' Vincenty Direct Calculation of Distance and Direction
#'
#' \code{distance} estimates the distance given a starting & ending latitude
#' and longitude. \cr \cr For general information on Vincenty's formula, see
#' e.g., \url{http://en.wikipedia.org/wiki/Vincenty's_formulae}. It states: \cr
#' \emph{Vincenty's formulae are two related iterative methods used in geodesy
#' to calculate the distance between two points on the surface of an spheroid,
#' developed by Thaddeus Vincenty in 1975. They are based on the assumption
#' that the figure of the Earth is an oblate spheroid, and hence are more
#' accurate than methods such as great-circle distance which assume a spherical
#' Earth.} \cr \cr \bold{Note:} this method assumes a locations are lat & lon
#' given in WGS 84.\cr\cr Direction, if requested, is the the initial bearing
#' (sometimes referred to as forward azimuth) for which one would follow as a
#' straight line along a great-circle arc from start to finish.\cr \cr
#' \bold{Note:} this will fail if there are NA's in the data.
#'
#'
#' @param lat1 a single value or vector of values representing latitude in
#' decimal degrees from -90 to 90 degrees. Alternatively, a data.frame or
#' matrix can be used here with each column representing lat1, lon1, lat2, lon2
#' (in that order).
#' @param lon1 a single value or vector of values representing longitude in
#' decimal degrees from -180 to 180 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param lat2 a single value or vector of values representing latitude in
#' decimal degrees from -90 to 90 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param lon2 a single value or vector of values representing longitude in
#' decimal degrees from -180 to 180 degrees. If NULL, lat1 is assumed to be a
#' matrix or data.frame.
#' @param bearing boolean value as to calculate the direction as well as the
#' distance.
#' @return Returns a data.frame with: \item{lon1}{the original longitude}
#' \item{lat1}{the original latitude} \item{lon2}{the destination longitude}
#' \item{lat2}{the destination latitude} \item{distance}{the distance used}
#' \item{bearing}{if requested, the bearing between the two points}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{destination}}
#' @references Vincenty, T. 1975. Direct and Inverse Solutions of Geodesics on
#' the Ellipsoid with application of Nested Equations. Survey Review, vol XXII
#' no 176. \url{http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf}
#' @source The source code for the distance algorithm here was modified from
#' \url{http://www.movable-type.co.uk/scripts/latlong-vincenty.html}.\cr \cr
#' Distances were validated against Geoscience Australia calculations
#' (\url{http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp}).\cr \cr
#' Bearings were from multiple sources including
#' \url{http://williams.best.vwh.net/avform.htm#Crs}.
#' @examples
#'
#'
#' #get the distance of 1 degree longitude at each 5 degrees latitude from -90 to 90
#' distance(lat1=seq(-90,90,5),lon1=rep(0,37),lat2=seq(-90,90,5),lon2=rep(1,37),bearing=TRUE)
#'
#'
#' @export
#' @useDynLib SDMTools Dist
# This is a fix for SDM tool distance function due to floating point operation.
EC_Distance <- function(lat1, lon1=NULL, lat2=NULL, lon2=NULL, bearing=FALSE) {
    if (is.data.frame(lat1) | is.matrix(lat1)) { #if input is matrix or data.frame... break it out to individual vectors
        lat1 = as.matrix(lat1); if (ncol(lat1)!=4) stop('incorrect lat/lon inputs... must be matrix with 4 columns or 4 vectors')
        lon2=lat1[,4]; lat2=lat1[,3]; lon1=lat1[,2]; lat1=lat1[,1] #break out individual columns
    } else if (!is.null(lat2) & !is.null(lon1) & !is.null(lon2)) {
        if (!all(c(length(lat2),length(lon1),length(lon2))==length(lat1))) stop('inputs must all be of same length')
    } else { stop('inappropriate inputs... see helpfile') }
    if (any(c(lon1,lon2) < -180.0001) | any(c(lon1,lon2) > 180.0001)) stop('lon must be decimal degrees between -180 & 180')
    if (any(c(lat1,lat2) < -90.0001) | any(c(lat1,lat2) > 90.0001)) stop('lat must be decimal degrees between -90 & 90')
    #cycle through and output the new data
    out = data.frame(lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2)
    out$distance = round(.Call('Dist',out$lat1,out$lon1,out$lat2,out$lon2,PACKAGE='SDMTools'),2) #round to the nearest mm
    if (bearing) { #if requested, calculate bearing
        lat1=lat1*pi/180;lat2=lat2*pi/180;lon1=lon1*pi/180;lon2=lon2*pi/180 #convert to radians
        brng = atan2(sin(lon2-lon1)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon1-lon2)) #estimate bearing
        out$bearing = ((brng*180/pi)+360)%%360 #convert to bearing in degrees
    }
    #return the output
    return(out)
}


########################################################
# Patch Biomod2 sample.factor.levels function
#  taken from: biomod2-3.3-7 : https://github.com/cran/biomod2/blob/master/R/sample.factor.levels.R#L70
########################################################

EC_SampleFactorLevels <- function(x, mask.out = NULL, mask.in = NULL){
  ## make some checking of given parameters
  ## TODO(damien)
  if(inherits(x, 'Raster')){
    fact.level.cells <- EC_SampleFactorLevelsRaster(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else if(inherits(x, 'data.frame')){
    fact.level.cells <- EC_SampleFactorLevelsDataFrame(x, mask.out = mask.out, mask.in = mask.in)
    return(fact.level.cells)
  } else {
    warning(paste0("\nunsupported input data.",
                   "\nx should be a Raster* object or a data.frame.",
                   "\n NULL returned"))
    return(NULL)
  }
}

EC_SampleFactorLevelsRaster <- function (x, mask.out = NULL, mask.in = NULL) 
{
  fact.var <- which(is.factor(x))
  if (any(fact.var)) {
    fact.level.cells <- sapply(fact.var, function(f) {
      selected.cells <- NULL
      fact.level.original <- unlist(raster::levels(subset(x, 
                                                          f)))
      fact.level <- fact.level.original
      cat("\n> fact.level for", names(x)[f], ":\t", paste(fact.level, 
                                                          names(fact.level), sep = ":", collapse = "\t"))
      if (!is.null(mask.out)) {
        fact.levels.sampled <- unlist(levels(as.factor(mask(subset(x, 
                                                                   f), mask.out))))
        attr(fact.levels.sampled, "names") <- attr(fact.level.original, 
                                                   "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, 
            "have already been sampled")
        fact.level <- fact.level[!is.element(fact.level,
                                             fact.levels.sampled)]
      }
      if (length(fact.level)) {
        if (!is.null(mask.in)) {
          for (mask.in.id in 1:length(mask.in)) {
            if (length(fact.level)) {
              x.f.masked <- as.factor(mask(subset(x, 
                                                  f), mask.in[[mask.in.id]]))
              x.f.levels <- unlist(levels(x.f.masked))
              attr(x.f.levels, "names") <- attr(fact.level.original, 
                                                "names")[x.f.levels]
              fact.levels.in.m.in <- fact.level[is.element(fact.level, 
                                                           x.f.levels)]
              if (length(fact.levels.in.m.in)) {
                cat("\n - levels", fact.levels.in.m.in, 
                    "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, 
                                                           function(fl) {
                                                             sample(which(x.f.masked[] == fl), 
                                                                    1)
                                                           }))
                fact.level <- fact.level[!is.element(fact.level, 
                                                     fact.levels.in.m.in)]
              }
            }
          }
        }
        if (length(fact.level)) {
          cat("\n - levels", fact.level, "will be sampled in the original raster")
          selected.cells <- c(selected.cells, sapply(fact.level, 
                                                     function(fl) {
                                                       sample(which(subset(x, f)[] == fl), 1)
                                                     }))
        }
      }
      return(selected.cells)
    })
    fact.level.cells <- as.numeric(fact.level.cells[-which(sapply(fact.level.cells, is.null))])
    return(fact.level.cells)
  }
  else {
    return(NULL)
  }
}

EC_SampleFactorLevelsDataFrame <- function(x, mask.out = NULL, mask.in = NULL){
  ## identificate the factorial variables
  fact.var <- which(sapply(x, is.factor))
  ## check if some factorial variables are in the input data
  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      ## initialize the list of cells that are selected
      selected.cells <- NULL
      ## get the levels of the factor on the full dataset
      fact.level.original <- levels(x[, f])
      fact.level <- fact.level.original
      cat("\n> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), fact.level, sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the fector that have been already sampled
        fact.levels.sampled <- unique(na.omit(as.character(x[mask.out[, 1], f])))
        ## remove already sampled points from candidates
        x[mask.out[, 1], ] <- NA
        #         ## update levels names (lost during mask conversion)
        #         attr(fact.levels.sampled, "names") <- attr(fact.level.original, "names")[fact.levels.sampled]
        cat("\n - according to mask.out levels", fact.levels.sampled, "have already been sampled")
        ## update the list of factor levels to sample
        fact.level <- setdiff(fact.level, fact.levels.sampled)
      }
      if(length(fact.level)){
        ## try first to sample factors in the given masks
        if(!is.null(mask.in)){ ## list of mask we want to sample in (order matter!)
          for(mask.in.id in 1:ncol(mask.in)){
            ## check that some levels remains to be sampled
            if(length(fact.level)){
              ## update the masked version of the factorial raster
              x.f.masked <- as.character(x[, f])
              x.f.masked[!mask.in[, mask.in.id]] <- NA
              x.f.levels <- unique(na.omit(x.f.masked))
              ## get the list of levels that coulb be sampled in this mask
              fact.levels.in.m.in <- intersect(fact.level, x.f.levels)
              if(length(fact.levels.in.m.in)){
                cat("\n - levels", fact.levels.in.m.in, "will be sampled in mask.out", mask.in.id)
                selected.cells <- c(selected.cells, sapply(fact.levels.in.m.in, function(fl){
                  candidate.cells <- na.omit(which(x.f.masked[] == fl))
                  selected.cell <- NULL
                  if(length(candidate.cells) == 1){ ## single candiate cell
                    selected.cell <- candidate.cells
                  } else { ## multi candidate cells
                    selected.cell <- sample(candidate.cells, 1)
                  }
                  return(selected.cell)
                }))
                ## update the list of factor levels to sample
                fact.level <- setdiff(fact.level, fact.levels.in.m.in)
              }
            } 
          } ## end loop over mask.in
        }
        ## @note if some levels remains unsampled then we will take a random value of
        ## them in the full dataset => !! this should be tricky if mask.in arg is given
        ## because the value will be picked out of mask.in but is necessary to 
        ## ensure that models will run smoothly
        if(length(fact.level)){
          cat("\n - levels", fact.level, "will be sampled in the original data.frame")
          selected.cells <- c(selected.cells, sapply(fact.level, function(fl){
            candidate.cells <- na.omit(which(x[, f] == fl))
            selected.cell <- NULL
            if(length(candidate.cells) == 1){ ## single candiate cell
              selected.cell <- candidate.cells
            } else { ## multi candidate cells
              selected.cell <- sample(candidate.cells, 1)
            }
            return(selected.cell)
          }))
        }
      }
      return(selected.cells)
    })))
    return(unique(fact.level.cells))
  } else { ## no factorial variable
    return(NULL)
  }
}

#assignInNamespace("SampleFactorLevels",SampleFactorLevels, ns="biomod2")
#assignInNamespace("EC_GridInfo",EC_GridInfo, ns="SDMTools")
#assignInNamespace("distance",EC_Distance, ns="SDMTools")

writeLines('Done')