# SignatureDataAnalysis.R

#Source associated scripts for data
source("./scripts/DataGeneration.R")

#Check which genes are missing in the signatures DB
which(!(GLYCOLYSIS_GENES %in% rownames(signatures)))
which(!(GLYCOSYLATION %in% rownames(signatures)))
which(!(RIBOSOMAL_S %in% rownames(signatures)))
# There seems to be reasonable coverage here

####-------- Correlation Analysis -------####

# Compile all systems together for comparison
AllGenes <- c(GLYCOLYSIS_GENES, GLYCOSYLATION, RIBOSOMAL_S)
subset <- AllGenes[AllGenes %in% rownames(signatures)]
systems <- c(rep("Glycolysis", length(GLYCOLYSIS_GENES)),
             rep("Glycosylation", length(GLYCOSYLATION)),
             rep("Ribosome", length(RIBOSOMAL_S)))
systems <- data.frame(system = systems, row.names = AllGenes,
                      stringsAsFactors = F)
# Select the rows featuring our investigated genes
Interest <- as.matrix(signatures[subset,])
# Empty data frame with genes to be coerced to rownames
GeneCorrelation <- data.frame(gene = subset, stringsAsFactors = F)
for(i in 1:length(subset)){
  # Correlation between row i and every other, including i which is 1
  corrs <- apply(Interest, 1L, cor, Interest[i,])
  # Add this to column with i gene name
  GeneCorrelation[,subset[i]] <- corrs
  # Set the diagonal to NA, since it is a gene and itself
  GeneCorrelation[i,i+1] <- NA
}

# New column with labels the system the gene belongs
GeneCorrelation$System <- systems[GeneCorrelation$gene,1]

# Visually check the mean and sd corr by gene system
corMeans <- apply(GeneCorrelation[,c(-1,-ncol(GeneCorrelation))], 2L, 
                  tapply, GeneCorrelation$System, mean, na.rm=T)
corSD <- apply(GeneCorrelation[,c(-1,-ncol(GeneCorrelation))], 2L, 
               tapply, GeneCorrelation$System, sd, na.rm=T)

# Create a df with just the corr between pairs of genes rowise
GeneByGene <- gather(GeneCorrelation, key = gene2, value = corr,
                     2:(ncol(GeneCorrelation)-1), na.rm = T) 
# Find the unique rows since there will be duplicates of every pair
uniqueRows <- lapply(strsplit(paste(GeneByGene$gene, GeneByGene$gene2), split=" "), 
                     sort)
GeneByGene <- GeneByGene[which(!duplicated(uniqueRows)),]
# Add system for gene 2
GeneByGene$System2 <- GeneCorrelation[GeneByGene$gene2, "System"]
# Create a match stat which will be the basis for training
GeneByGene$Match <- ifelse(GeneByGene$System == GeneByGene$System2,
                           1, 0)

set.seed(pi)
# Select rows for training and testing
trainRows <- sample(1:nrow(GeneByGene), 2000)
train <- GeneByGene[trainRows,]
test <- GeneByGene[-trainRows,]
# Create an lm based on a binomial model
model <- glm(Match ~ corr, family = binomial(link = "logit"), train)
summary(model)
# We can see that this is apparently significant, though the dataset is large

# Compare to our testing set values
fittedValues <- predict(model, newdata=select(test, corr), type='response')
fittedValues <- ifelse(fittedValues > 0.5, 1, 0)
errors <- mean(fittedValues != test$Match)
FP <- mean(fittedValues == 1 & test$Match == 0)
# 0.0608
FN <- mean(fittedValues == 0 & test$Match == 1)
# 0.161...
# We see that this initial model has an accuracy of only 77.8%


####-------- Neural Network Three values --------- ####

geneMatch <- data.frame()
for (i in 1:(nrow(Interest) - 2)){
  # Search down from the current index
  geneMatrix <- Interest[-(1:i),]
  # Transpose to compare a row to columns
  geneMatrix <- t(geneMatrix)
  geneInterest <- Interest[i,]
  # Mismatches are set to -1
  geneMatrix[geneMatrix != geneInterest] <- -1
  # Turn this back into a df and add the gene names as columns
  geneMatrix <- as.data.frame(t(geneMatrix))
  geneMatrix$gene1 <- rep(rownames(Interest)[i], nrow(geneMatrix))
  geneMatrix$gene2 <- rownames(Interest)[-(1:i)]
  # Bind this to the master list
  geneMatch <- rbind(geneMatch, geneMatrix)
}
# Last row of the df fails add it seperately
i <- i +1
geneMatrix <- Interest[i,]
geneMatrix[geneMatrix != Interest[i+1,]] <- -1
geneMatrix <- as.data.frame(t(geneMatrix))
geneMatrix$gene1 <- rownames(Interest)[i]
geneMatrix$gene2 <- rownames(Interest)[i+1]
geneMatch <- rbind(geneMatch, geneMatrix)

# Add the system based on the gene name and check whether they match
geneMatch$system1 <- systems[geneMatch$gene1, 1]
geneMatch$system2 <- systems[geneMatch$gene2, 1]
geneMatch$Match <- ifelse(geneMatch$system1 == geneMatch$system2, 1, 0)

set.seed(pi)
# Select the rows for training and testing
trainRows <- sample(1:nrow(geneMatch), 2000)
train <- geneMatch[trainRows,] %>%
  select(-gene1, -gene2, -system1, -system2)
test <- geneMatch[-trainRows,] %>%
  select(-gene1, -gene2, -system1, -system2, -Match)
matches <- geneMatch[-trainRows, "Match"]

# Create the formula which includes every column except the match
formula <- as.formula(paste("Match ~ ", 
                            paste(colnames(train)[-ncol(train)],
                                  collapse = "+"), sep = ""))
# Create a neural network, with two hidden layers
model <- neuralnet(formula = formula, train, hidden = 2)
# Check the values predicted based on this new model
predictions <- compute(model, test)
predictions <- ifelse(predictions$net.result > 0.5, 1, 0)

errors <- mean(predictions != geneMatch$Match[-trainRows])
FP <- mean(predictions == 1 & geneMatch$Match[-trainRows] == 0)

##### ------- MCA For Dimensionality Reduction ------- #####
binaryMatrix <- apply(as.data.frame(Interest), 2, as.factor)
geneMCA <- MCA(binaryMatrix, ncp = 1, graph = F)
Projection <- geneMCA$ind$coord
names <- names(geneMCA$ind$contrib)
result <- Ckmeans.1d.dp(Projection, k=c(1,7))
cluster <- result$cluster
compare <- data.frame(gene1 = rep(names, length(names)),
                      gene2 = as.vector(sapply(names, rep, list(length(names)))),
                      cluster1 = rep(cluster, length(names)),
                      cluster2 = as.vector(sapply(cluster, rep, list(length(names)))),
                      stringsAsFactors = F) %>%
  mutate(gMatch = cluster1 == cluster2,
         system1 = systems[gene1, 1],
         system2 = systems[gene2, 1],
         rMatch = system1 == system2) %>%
  filter(gene1 < gene2)
mean(compare$rMatch != compare$gMatch)
sum(compare$rMatch == 0 & compare$gMatch == 1) / sum(compare$rMatch == 0)

# [END]