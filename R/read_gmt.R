# read_gmt.R

gmtToDF <- function(filename) {
  #
  # Reads in from filename and returns a DF with column 1 as gene names, other
  # columns logical as to whether that gene was successful in the experiment
  
  genes <- read_file(filename)
  genes <- strsplit(genes, "\n")
  genes <- lapply(genes, strsplit, split = "\t")[[1]]
  
  ExpName <- sapply(genes, function(x) x[[1]]) # Get standard name for experiment
  # Remove the experiment name and comment entries
  genes <- sapply(genes, function(x) x[-1:-2])

  # Creates vector of groups based on length of each experiment group
  groups <- unlist(mapply(function(x,y) rep(x, length(y)),
                    1:length(genes),
                    genes))
  # Adding unique to df creation prevents an error in spread from error columns
  # with duplicate genes from a single dataset
  GeneMembership <- unique(data.frame(gene = unlist(genes),
                        ExpID = groups,
                        stringsAsFactors = F)) %>%
    mutate(True = T) %>% # Sets a column to true, the case for each experiment
    spread(ExpID, True, fill = F) # Turns genes into rows, with logical column for each experiment
  
  colnames(GeneMembership) <- c("gene", ExpName) # Set as column names
  
  return(GeneMembership)
}

# [END]