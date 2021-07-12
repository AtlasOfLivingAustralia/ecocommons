#' Computes the Area/Extent of Occuppancy (AOO/EOO) for species of interest
#'
#' @param occur
#' @param species
#'
#' @return
#' @export
#'
#' @examples

EC_AreaExtentOccupancy_test <- function(occur, species) {
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
