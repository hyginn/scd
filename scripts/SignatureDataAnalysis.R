# SignatureDataAnalysis.R

# Create with gene symbols for the systems of interest
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
rownames(signatures) <- signatures[,1]
signatures <- signatures[,-1]
signatures[is.na(signatures)] <- F

which(!(GLYCOLYSIS_GENES %in% rownames(signatures)))
which(!(GLYCOSYLATION %in% rownames(signatures)))
which(!(RIBOSOMAL_S %in% rownames(signatures)))

AllGenes <- c(GLYCOLYSIS_GENES, GLYCOSYLATION, RIBOSOMAL_S)
for(i in 1:(length(AllGenes)-1)){
  for(j in i:length(AllGenes)){
    
  }
}