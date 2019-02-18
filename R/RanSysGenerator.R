# RanSysGenerator.R


library(ggplot2)

set.seed(473)

RanSysGenerator <- function(empd, MINSIZE = 10, MAXSIZE = 30, 
                            nsys = 100, nfeatures = 5000){
  #
  # Use an empirical distribution of probabilities to generate a binary matrix
  # empd: empirical distribution of probabilities to draw from
  # MINSIZE: the minimum number of "genes" in a system
  # MAXSIZE: maximum number of genes present in any given system
  # nsys: number of systems required to be made
  # nfeatures: number of synthetic features required to be generated
  if (MINSIZE >= MAXSIZE) {
    stop(sprintf("PANIC: MINSIZE (%d) is not less then MAXSIZE (%d)",
                 MINSIZE,
                 MAXSIZE))
  }
  # Select the new probabilities of a gene in a system having any given feature
  # Each value represents a single system-feature pair
  probs <- sample(empd, nsys*nfeatures, replace = T)
  # Randomly determine the sizes of our genes
  SysSizes <- sample(MINSIZE:MAXSIZE, nsys, replace = T)
  
  # We will now use the probabilities to draw from a binom distribution
  # For a gene to contain a certain feature
  
  # Prepare data to use vectorized rbinom
  totalGenes <- matrix(ncol = nfeatures)
  for(i in 1:nsys){
    # Create a matix by system
    # Each system will repeat the probability of a certain feature for each gene
    # columns = features, rows = genes
    temp <- rep(probs[(nfeatures*(i-1) + 1):(nfeatures*i)], SysSizes[i])
    temp <- matrix(temp, nrow = SysSizes[i], byrow = T)
    totalGenes <- rbind(totalGenes, temp)
  }
  # First row is NA, remove
  totalGenes <- totalGenes[-1,]
  # Take binom of 1 for ech probability present in totalGenes
  Attributes <- rbinom(length(totalGenes), 1, totalGenes)
  # Transform into df
  Attributes <- as.data.frame(matrix(Attributes, ncol = nfeatures))
  # Add the system labels
  Attributes$System <- unlist(mapply(function (x,y) rep(x, y), 1:nsys, SysSizes))
  return(Attributes[sample(nrow(Attributes)), ]) # Return df with rows
                                                 # randomized
}

# [END]
