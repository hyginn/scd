# SignatureDataAnalysis.R

# Create vectors with gene symbols for the systems of interest
GLYCOLYSIS_GENES <- read.delim("./data/glycolysis.tab", stringsAsFactors = F)$Symbol
RIBOSOMAL_S <- read.delim("./data/ribosomal.tab", stringsAsFactors = F)$Approved.symbol
GLYCOSYLATION <- read.delim("./data/glycosylation.tab", stringsAsFactors = F)$ProteinHuman
GLYCOSYLATION <- toupper(GLYCOSYLATION)[GLYCOSYLATION != "?"]

# Read in all our data sets from MSigDB into dfs
Hallmark <- gmtToDF("./data/Hallmark.gmt")
Curated <- gmtToDF("./data/Curated.gmt")
CompSet <- gmtToDF("./data/CompGeneSet.gmt")

# Join these frames together
signatures <- full_join(Hallmark, Curated, by = "gene")
signatures <- full_join(signatures, CompSet, by = "gene")
# Add genes as rownames instead of a column
rownames(signatures) <- signatures[,1]
signatures <- signatures[,-1]
# Full join puts missing data to NA, replace with False instead
signatures[is.na(signatures)] <- F

#Check which genes are missing in the signatures DB
which(!(GLYCOLYSIS_GENES %in% rownames(signatures)))
which(!(GLYCOSYLATION %in% rownames(signatures)))
which(!(RIBOSOMAL_S %in% rownames(signatures)))
# There seems to be reasonable coverage here


# Compile all systems together for comparison
AllGenes <- c(GLYCOLYSIS_GENES, GLYCOSYLATION, RIBOSOMAL_S)
# Select the rows featuring our investigated genes
Interest <- as.matrix(signatures[AllGenes,])
# Empty data frame with genes to be coerced to rownames
GeneCorrelation <- data.frame(gene = AllGenes, stringsAsFactors = F)
for(i in 1:length(AllGenes)){
  # Correlation between row i and every other, including i which is 1
  corrs <- apply(Interest, 1L, cor, Interest[i,])
  # Add this to column with i gene name
  GeneCorrelation[,AllGenes[i]] <- corrs
  # Set the diagonal to NA, since it is a gene and itself
  GeneCorrelation[i,i+1] <- NA
}

# New column with labels the system the gene belongs
GeneCorrelation$System <- c(rep("Glycolysis", length(GLYCOLYSIS_GENES)),
                            rep("Glycosylation", length(GLYCOSYLATION)),
                            rep("Ribosome", length(RIBOSOMAL_S)))

# Visually check the mean and sd corr by gene system
corMeans <- apply(GeneCorrelation[,c(-1,-ncol(GeneCorrelation))], 2L, 
                  tapply, GeneCorrelation$System, mean, na.rm=T)
corSD <- apply(GeneCorrelation[,c(-1,-ncol(GeneCorrelation))], 2L, 
               tapply, GeneCorrelation$System, sd, na.rm=T)

# Create a df with just the corr between pairs of genes rowise
GeneByGene <- gather(GeneCorrelation, key = gene2, value = corr,
                     2:(ncol(GeneCorrelation)-1), na.rm = T) %>%
  unique()
# Add the system of the second gene to the df for comparison
rownames(GeneCorrelation) <- GeneCorrelation$gene
GeneByGene$System2 <- GeneCorrelation[GeneByGene$gene2, "System"]
# Create a match stat which will be the basis for training
GeneByGene$Match <- ifelse(GeneByGene$System == GeneByGene$System2,
                           1, 0)


set.seed(pi)
# Select rows for training and testing
trainRows <- sample(1:nrow(GeneByGene), 4000)
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
# 0.053768...
FN <- mean(fittedValues == 0 & test$Match == 1)
# 0.17224...
# We see that this initial model has an accuracy of only 77.4%

# [END]