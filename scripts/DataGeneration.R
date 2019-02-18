####-------- Data Generation -------####

# DataGeneration.R
source("./R/include.R")
source("./R/read_gmt.R")

# Create vectors with gene symbols for the systems of interest
GLYCOLYSIS_GENES <- read.delim("./data/glycolysis.tab", stringsAsFactors = F)$Symbol
RIBOSOMAL_S <- read.delim("./data/ribosomal.tab", stringsAsFactors = F)$Approved.symbol
GLYCOSYLATION <- read.delim("./data/glycosylation.tab", stringsAsFactors = F)$ProteinHuman
GLYCOSYLATION <- toupper(GLYCOSYLATION)[GLYCOSYLATION != "?"]
SONIC <- read.delim("./data/shh.tab", stringsAsFactors = F)$Gene
RNA_POL <- read.delim("./data/RNAPol.tab", stringsAsFactors = F)$Approved.symbol

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
