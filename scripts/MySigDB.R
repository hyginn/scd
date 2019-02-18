# MySigDB.R
library(readr)
library(dplyr)

read_gmt <- function(filename, institute = "broadinstitute") {
  #
  # institute: The source of the data, seperates experiments in gmt file
  #
  # Reads in from a file and returns a DF with column 1 as gene names, other
  # columns logical as to whether that gene was successful in the experiment
  
  genes <- read_file(filename)
  genes <- strsplit(genes, "[\t\n]")
  
  # Finds the start index of each of the experiments, using source website for search
  indexes <- which(sapply(genes, gregexpr, pattern = institute) > 0)
  # Sorts each index into a group based on which experiment it belongs
  groups <- sapply(1:length(genes[[1]]), function(x) sum(x > indexes))
  TotalDF <- data.frame(gene = genes[[1]],
                        ExpID = groups,
                        stringsAsFactors = F) 
  # Includes the experiment names and source website
  GeneMembership <-
    TotalDF[c(-indexes,-(indexes - 1)), ] %>% # Removes source website and experiment name
    mutate(True = T) %>% # Sets a column to true, the case for each experiment
    spread(ExpID, True, fill = F) # Turns genes into rows, with logical column for each experiment
  
  ExpName <- genes[[1]][indexes - 1] # Get standard name for experiment
  colnames(GeneMembership) <- c("gene", ExpName) # Set as column names
  
  return(GeneMembership)
}

