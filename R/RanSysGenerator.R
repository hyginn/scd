# RanSysGenerator.R


library(ggplot2)

set.seed(473)

RanSysGenerator <- function(means, sd, MINSIZE = 10, MAXSIZE = 30){
  #
  # MINSIZE: the minimum number of "genes" in a system
  #
  # Creates N systems with attributes centered around means with a uniform
  # standard deviation of sd
  if (MINSIZE >= MAXSIZE) {
    stop(sprintf("PANIC: MINSIZE (%d) is not less then MAXSIZE (%d)",
                 MINSIZE,
                 MAXSIZE))
  }
  Attributes <- data.frame() # Create the observation DF
  for(i in 1:length(means)){
    # Create a random number of observations following a system specific ND
    value <- rnorm(sample(MINSIZE:MAXSIZE, 1), means[i], sd)
    system <- rep(i, length(value)) # Create labels with the system for
                                    # validation
    Attributes <- rbind(Attributes, cbind(system, value)) # Store in Observation DF
  }

  return(Attributes[sample(nrow(Attributes)), ]) # Return df with rows
                                                 # randomized
}


GeneratedSys <- RanSysGenerator(c(3, 8, 20), 1) # Create our random data
clusters <- kmeans(GeneratedSys$value, 3) # Generic k-means clustering on values
# Label each row with the mean of its assigned cluster
GeneratedSys$clusterMean <- clusters$centers[clusters$cluster]

# Plot the points on a line, colouring them based on their assigned mean
ggplot(data = GeneratedSys, # Create a
       aes(x = value, y = 0, color = as.factor(cluster))) +
  geom_point(size = 3) + ylab("") +
  scale_color_discrete(name = "Cluster Mean")

# Find the number of unique cluster means for each of our determined systems
tapply(GeneratedSys$cluster, as.character(GeneratedSys$system),
       function(x) length(unique(x)))
#This should return all 1, if our clustering was done correctly

# [END]
