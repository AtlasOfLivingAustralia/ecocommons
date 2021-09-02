#' Use BIOMOD_FormatingData function from 'biomod2' package to format user
#' input data
#' 
#' It generates pseudo absence points if true absence data are not available or
#' adds pseudo absence data to an existing absence dataset
#' 
#' @param true.absen 
#' @param pseudo.absen.points 
#' @param pseudo.absen.strategy 
#' @param pseudo.absen.disk.min 
#' @param pseudo.absen.disk.max 
#' @param pseudo.absen.sre.quant 
#' @param climate.data 
#' @param occur 
#' @param species.name 
#' @param save.pseudo.absen 
#' @param save.env.absen 
#' @param save.env.occur 
#' @param generate.background.data 
#' @param species_algo_str 
#'

#' @export EC_format_biomod2
#' @importFrom biomod2 BIOMOD_FormatingData
#'

EC_format_biomod2 <- function (true.absen               = NULL,
                               pseudo.absen.points      = 0,
                               pseudo.absen.strategy    = 'random',
                               pseudo.absen.disk.min    = 0,
                               pseudo.absen.disk.max    = NULL,
                               pseudo.absen.sre.quant   = 0.025,
                               climate.data             = NULL,
                               occur                    = NULL,
                               species.name             = NULL,
                               save.pseudo.absen        = TRUE,
                               save.env.absen           = TRUE,
                               save.env.occur           = TRUE,
                               generate.background.data = FALSE,
                               species_algo_str         = NULL) {
  
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
    absen = data.frame(lon=numeric(0), lat=numeric(0))  # create an empty data frame for background points
    pseudo.absen.rep = 1  # generate pseudo-absence points
    if (!save.pseudo.absen) {
      pseudo.absen.rep = 0
    }
  } else {
    absen <- true.absen
    if (!is.null(climate.data) && nrow(true.absen) > 0) { # ensure true absence is in same projection as climate
      absen <- EC_data_projection(true.absen, climate.data)
    }
    
    pseudo.absen.rep = 0 # do not generate pseudo absence point when true absence points are available
    pseudo.absen.strategy = 'none'
    pseudo.absen.points = nrow(absen)
    cat("No pseudo absence point is generated")
    }
  
  # generate background data as pseudo absence points
  if (pseudo.absen.strategy != 'none' & generate.background.data) {
    biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))
    myBackgrdData <-
      biomod2::BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                                    expl.var  = climate.data,
                                    resp.name = species.name,
                                    PA.nb.rep = pseudo.absen.rep,
                                    PA.nb.absences = pseudo.absen.points,
                                    PA.strategy = pseudo.absen.strategy,
                                    PA.dist.min = pseudo.absen.disk.min,
                                    PA.dist.max = pseudo.absen.disk.max,
                                    PA.sre.quant = pseudo.absen.sre.quant)
    
    # get background data as absence data
    colnames(myBackgrdData@coord) <- c('lon', 'lat')
    absen <- myBackgrdData@coord[c(which(is.na(myBackgrdData@data.species))), c('lon', 'lat')]
    
    # do not generate pseudo absence in next call to BIOMOD_FormatingData
    pseudo.absen.rep = 0
    pseudo.absen.strategy = 'none'
  }
  
  biomod.data <- rbind(occur[,c("lon", "lat")], absen[,c("lon", "lat")])
  biomod.data.pa <- c(rep(1, nrow(occur)), rep(0, nrow(absen)))
  
  myBiomodData <-
    biomod2::BIOMOD_FormatingData(resp.var  = biomod.data.pa,
                                  expl.var  = climate.data,
                                  resp.xy   = biomod.data,
                                  resp.name = species.name,
                                  PA.nb.rep = pseudo.absen.rep,
                                  PA.nb.absences = pseudo.absen.points,
                                  PA.strategy = pseudo.absen.strategy,
                                  PA.dist.min = pseudo.absen.disk.min,
                                  PA.dist.max = pseudo.absen.disk.max,
                                  PA.sre.quant = pseudo.absen.sre.quant)
  
  # save the pseudo absence points generated to file
  pa_filename <- EC_outfile_name (filename="pseudo_absences",
                                  id_str=species_algo_str, ext="csv")
  absenv_filename <- EC_outfile_name (filename="absence_environmental",
                                      id_str=species_algo_str, ext="csv")
  occenv_filename <- EC_outfile_name (filename="occurrence_environmental",
                                      id_str=species_algo_str, ext="csv")
  if (pseudo.absen.rep != 0) {
    pseudoAbsen <- myBiomodData@coord[c(which(is.na(myBiomodData@data.species))), c('lon', 'lat')]
    if (save.pseudo.absen & nrow(pseudoAbsen) > 0) {
      EC_write_csv (pseudoAbsen, pa_filename, rownames = FALSE)
    }
    
    # save the pseudo absence points with environmental variables
    if (save.env.absen) {
      EC_merge_save (climate.data, pseudoAbsen, species.name, absenv_filename)
    }
  }
  else if (nrow(absen) > 0) {
    # save true-absence/background data generated
    if (!is.null(true.absen)) {
      pa_filename = EC_outfile_name(filename="absence", id_str=species_algo_str,
                                    ext="csv")  # rename true-absence file
    }
    EC_write_csv(absen, pa_filename, rownames = FALSE)
    
    if (save.env.absen) {
      # save the true absence/background points with environmental data
      EC_merge_save (climate.data, absen, species.name, absenv_filename)
    }
  }
  
  if (save.env.occur) {
    # save the occurrence datasets with environmental variables
    EC_merge_save (climate.data, occur, species.name, occenv_filename)
  }
  
  return(myBiomodData)
}




#===================================================================
## Subfunctions, listed in the order in which they appear within EC_format_biomod2()


## This is a common function to merge and save the new generated data
## in a csv file

EC_merge_save <- function(env, csvdata, spname, ofname) {
  data <- cbind(csvdata, species = spname, extract(env, csvdata[c('lon','lat')]))
  
  EC_write_csv(data, ofname, rownames=FALSE)
}


