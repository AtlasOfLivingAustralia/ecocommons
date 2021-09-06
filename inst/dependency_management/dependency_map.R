# look at dependencies

# install.packages("miniCRAN")
library(miniCRAN)

# # install_github("AtlasOfLivingAustralia/ecocommons")
# library(ecocommons)
# ecocommons_deps <- pkgDep("ecocommons")  # fails

# # take code from depgraph to make this work
# remotes::install_github("crsh/depgraph")
# library("depgraph")
# plot_dependency_graph("ecocommons")

pkg <- devtools::as.package("ecocommons")
dependencies <- unlist(strsplit(pkg$imports, split = "\n"))
dependencies <- gsub("\\n| \\(.+\\)|,", "", dependencies)
dependency_list <- lapply(dependencies, pkgDep)
names(dependency_list) <- dependencies

# remove self names from dependency_list
dependency_list2 <- lapply(dependencies, function(a){
  result <- dependency_list[[a]][dependency_list[[a]] != a]
  if(length(result) > 0){
    data.frame(from = a, to = result)
  }else{
    NULL
  }
})

# convert to df
# NOTE this df includes all sub-dependencies
dependency_df <- rbind(
  data.frame(from = "ecocommons", to = dependencies),
  do.call(rbind, dependency_list2))

# work out how many times each package is called
dependency_df <- merge(dependency_df, 
  as.data.frame(xtabs(~ to, data = dependency_df)))

# could add download size here if needed


# map only the ecocommons direct dependencies, get summary stats
# convert to wide format
all_dependencies <- unique(do.call(c, dependency_df))
dep_matrix <- do.call(rbind, 
  lapply(dependency_list, function(a){
    as.numeric(all_dependencies %in% a)
  }))
rownames(dep_matrix) <- dependencies
colnames(dep_matrix) <- all_dependencies

dep_stats <- data.frame(
  name = dependencies,
  n_packages = apply(dep_matrix, 1, sum))
rownames(dep_stats) <- NULL

col_sums <- apply(dep_matrix, 2, sum)
dep_stats$unique_packages <- unlist(lapply(
  split(dep_matrix, seq_len(nrow(dep_matrix))),
  function(a){sum(as.numeric(col_sums[a > 0] < 2))}))
dep_stats$prop_unique <- dep_stats$unique_packages / dep_stats$n_packages


# next step is to work out which functions are called most, and from which packages
# parse the information in @importFrom in roxygen2 text
all_functions <- lapply(
  list.files("./ecocommons/R"),
  function(a){
    text <- scan(
      file = paste0("./ecocommons/R/", a), 
      what = "character",
      sep = "\n")
    keep_lines <- which(grepl("#' @importFrom", text))
    if(length(keep_lines) > 0){
      raw_text <- gsub("#' @importFrom ", "", text[keep_lines])
      text_df <- do.call(rbind, lapply(
        strsplit(raw_text, " "),
        function(x){data.frame(package = x[1], function_name = x[c(2:length(x))])}
      ))
      text_df$EC_function <- gsub(".R$", "", a)
      return(text_df)
    }else{
      NULL
    }
  })
all_functions_df <- do.call(rbind, all_functions)
all_functions_df <- all_functions_df[
  order(all_functions_df$package, all_functions_df$function_name), ]
write.csv(all_functions_df, "ecocommons_imports.csv", row.names = FALSE)

n_functions <- apply(
  xtabs(~ function_name + package, data = all_functions_df),
  2, function(a){length(which(a > 0))})

dep_stats <- merge(dep_stats, 
  data.frame(name = names(n_functions), n_functions = n_functions),
  all.x = TRUE, all.y = FALSE)

write.csv(dep_stats, "dependency_stats.csv", row.names = FALSE)


## BELOW HERE AN ALTERNATIVE WORKFLOW - NOT REQUIRED BUT POTENTIALLY USEFUL
dependency_graph <- miniCRAN::makeDepGraph(
  pkg = dependencies, 
  suggests = TRUE)
  
# plot interactively
library(networkD3)
library(igraph)
dependency_df <- as.data.frame(get.edgelist(dependency_graph))
colnames(dependency_df) <- c("to", "from")

# alternative pkg-level information (non-functional for our use-case)
install.packages("pkggraph")
library(pkggraph)
pkggraph::init(local = TRUE)
test <- get_neighborhood("ecocommons") # fails

# using DependenciesGraphs - better, but only at the function level
devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)
library(ecocommons)
dep <- envirDependencies("package:ecocommons")
plot(dep)