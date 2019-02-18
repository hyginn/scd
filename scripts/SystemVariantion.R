# SystemVariation.R

# Attempt to investigate the variation found within our characterized systems 
# and develop a method to simulate this for 10s-100s of systems

#Source associated scripts for data
source("./scripts/DataGeneration.R")

# Compile all systems together for comparison
AllGenes <- c(GLYCOLYSIS_GENES, GLYCOSYLATION, RIBOSOMAL_S, SONIC, RNA_POL)
subset <- AllGenes[AllGenes %in% rownames(signatures)]
systems <- c(rep("Glycolysis", length(GLYCOLYSIS_GENES)),
             rep("Glycosylation", length(GLYCOSYLATION)),
             rep("Ribosome", length(RIBOSOMAL_S)),
             rep("SHH", length(SONIC)),
             rep("RNA Pol", length(RNA_POL)))
# Select the rows featuring our investigated genes
Interest <- signatures[subset,]
# Indentify their systems
Interest$system <- systems[AllGenes %in% rownames(signatures)]

# Find the probability of a gene being present in every dataset
# Depending on its system
probs <- apply(Interest[,1:(ncol(Interest)-1)],2, 
               function(x) tapply(x, Interest$system, mean))

summary(t(probs))
library(ggplot2)
ggplot()+ geom_density(aes(x = probs[1,]), color = "blue") +
  geom_density(aes(x = probs[2,]), color = "red") +
  geom_density(aes(x = probs[3,]), color = "yellow")
