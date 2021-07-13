#' Patch Biomod2 sample.factor.levels function
#' Need to check if it is the same; if yes, delete
#' Three functions: 'EC_SampleFactorLevels', 'EC_SampleFactorLevelsRaster'
#' and 'EC_SampleFactorLevelsDataFrame'
#'
#' @param x 
#' @param mask.out 
#' @param mask.in 
#'
#' 

EC_SampleFactorLevels <- function(x, mask.out = NULL, mask.in = NULL){
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
  fact.var <- which(sapply(x, is.factor))    # identificate the factorial variables

  if(any(fact.var)){ ## some factorial variables present
    fact.level.cells <- as.numeric(unlist(sapply(fact.var, function(f){
      selected.cells <- NULL  # initialize the list of cells that are selected
      fact.level.original <- levels(x[, f])
      fact.level <- fact.level.original
      cat("\n> fact.level for",  colnames(x)[f], ":\t", paste(1:length(fact.level), 
                                                              fact.level, sep = ":", collapse = "\t"))
      if(!is.null(mask.out)){ ## mask containing points that have already been sampled
        ## check the levels of the factor that have been already sampled
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
