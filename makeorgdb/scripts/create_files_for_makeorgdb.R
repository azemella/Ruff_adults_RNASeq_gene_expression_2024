# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyr")
library("dplyr")

# Read the NCBI transcript dataset
# Link NCBI gene identifiers to NCBI gene descriptions
cpu.anno <- read.csv(file.path(dir, "metadata/ncbi-transcript-dataset.csv"))
cSym <- cpu.anno[,c(1,2,10)]
colnames(cSym)[1] <- "GID"
colnames(cSym)[2] <- "SYMBOL"
colnames(cSym)[3] <- "GENENAME"
cSym <- cSym[!duplicated(cSym), ]
write.csv(as.data.frame(cSym), 
          file = "makeorgdb/files/csym.csv" , row.names = FALSE)

# Link NCBI gene identifiers to Ensembl identifiers
cEns <- cpu.anno[,c(1,17,18,19)]
colnames(cEns)[1] <- "GID"
colnames(cEns)[2] <- "ENSEMBL"
colnames(cEns)[3] <- "ENSEMBLTRANS"
colnames(cEns)[4] <- "ENSEMBLPROT"
cEns <- cEns[!duplicated(cEns), ]
write.csv(as.data.frame(cEns), 
          file = "makeorgdb/files/cens.csv" , row.names = FALSE)

# Read the NCBI gene dataset
# Link NCBI gene identifiers to gene types
cpu.anno2 <- read.csv(file.path(dir, "metadata/ncbi-gene-dataset.csv"))
cType <- cpu.anno2[,c(1,5)]
colnames(cType)[1] <- "GID"
colnames(cType)[2] <- "GENETYPE"
cType <- cType[!duplicated(cType), ]
write.csv(as.data.frame(cType), 
          file = "makeorgdb/files/ctype.csv" , row.names = FALSE)

# Link NCBI gene identifiers to NCBI transcript and protein identifiers
cTx <- cpu.anno[,c(1,3,5)]
colnames(cTx)[1] <- "GID"
colnames(cTx)[2] <- "TXID"
colnames(cTx)[3] <- "PROTID"
cTx <- cTx[!duplicated(cTx), ]
write.csv(as.data.frame(cTx), 
          file = "makeorgdb/files/ctx.csv" , row.names = FALSE)

# Read the EggNOG dataset and extract protein IDs and GO terms
# Link NCBI gene identifiers to NCBI protein identifiers and eggNOG-mapper GO annotations
cpu.go <- read.csv(file.path(dir, "eggNOG-mapper/eggNOG_functional_annotations/out.emapper.annotations.csv"))
cGo <- cpu.go[,c(1,10)]
cGo <- cGo[!grepl("-", cGo$GOs),]
cGo <- cGo %>% dplyr::mutate(GOs = strsplit(GOs, ',')) %>% tidyr::unnest(cols = c(GOs))
cGo <- cGo[!duplicated(cGo), ]
colnames(cGo)[1] <- "PROTID"
colnames(cGo)[2] <- "GO"
cName <- cpu.anno[,c(1,5)]
colnames(cName)[1] <- "GID"
colnames(cName)[2] <- "PROTID"
cGOs <- merge(cName, cGo, by = "PROTID")
cGOs$EVIDENCE <- "ISO"
cGOs <- cGOs[-1]
cGOs <- cGOs[!duplicated(cGOs), ]
write.csv(as.data.frame(cGOs), 
          file = "makeorgdb/files/cgo.csv" , row.names = FALSE)

# Link NCBI gene identifiers to NCBI protein identifiers and eggNOG-mapper KEGG annotations
cKegg <- cpu.go[,c(1,12)]
cKegg <- cKegg[!grepl("-", cKegg$KEGG_ko),]
cKegg <- cKegg %>% dplyr::mutate(KEGG_ko = strsplit(KEGG_ko, ',')) %>% tidyr::unnest(cols = c(KEGG_ko))
cKegg <- cKegg[!duplicated(cKegg), ]
colnames(cKegg)[1] <- "PROTID"
colnames(cKegg)[2] <- "KEGG"
cKegg <- merge(cName, cKegg, by = "PROTID")
cKegg <- cKegg[-1]
cKegg <- cKegg[!duplicated(cKegg), ]
write.csv(as.data.frame(cKegg), 
          file = "makeorgdb/files/ckegg.csv", row.names = FALSE)
