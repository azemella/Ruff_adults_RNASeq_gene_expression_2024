# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("org.Cpugnax.eg.db")
library("clusterProfiler")

### Run for each Morph-Tissue comparison (36 comparisons in total)
### Modify lines #17, #31 and #97 accordingly

# Read the significant results file for the comparison
genes_ora <- read.csv(file.path(dir, 
                     "DESeq2_gene_expression_analysis/DESeq2_results/tissue/significant_results/DGE_DESeq2_satellite_vs_independent_Blo_SR.csv"))

# Extract up-regulated genes
upGenes <- genes_ora$Gene_symbol[ genes_ora$padj < 0.05 & 
                                  !is.na(genes_ora$padj) &
                                  genes_ora$log2FoldChange > 0 ]

# Extract down-regulated genes
downGenes <- genes_ora$Gene_symbol[ genes_ora$padj < 0.05 & 
                                    !is.na(genes_ora$padj) &
                                    genes_ora$log2FoldChange < 0 ]

# Read the universe background file
universe_background <- read.csv(file.path(dir, 
                      "DESeq2_gene_expression_analysis/DESeq2_results/tissue/full_results/DGE_DESeq2_satellite_vs_independent_Blo_FR.csv"))

# Extract universe genes
universe_genes <- universe_background[, 1]
universe_genes <- as.character(universe_genes)

# Perform Gene Ontology (GO) enrichment analysis for up-regulated genes
upgenes_ora <- enrichGO(gene = upGenes,
                        OrgDb = org.Cpugnax.eg.db,
                        keyType = "SYMBOL",
                        minGSSize = 10,
                        maxGSSize = 750,
                        ont = "ALL",
                        universe = universe_genes,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = TRUE)

# Perform Gene Ontology (GO) enrichment analysis for down-regulated genes
downgenes_ora <- enrichGO(gene = downGenes,
                          OrgDb = org.Cpugnax.eg.db,
                          keyType = "SYMBOL",
                          minGSSize = 10,
                          maxGSSize = 750,
                          ont = "ALL",
                          universe = universe_genes,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          readable = TRUE)

# Convert results to data frames
results_up_ora <- as.data.frame(upgenes_ora)
results_down_ora <- as.data.frame(downgenes_ora)

# Add column for up-regulated and down-regulated genes
results_up_ora$GeneList <- 'Up-regulated'
results_down_ora$GeneList <- 'Down-regulated'

# Combine results for up-regulated and down-regulated genes
results_ora_go <- rbind(results_up_ora, results_down_ora)

# Remove unwanted columns
results_ora_go <- results_ora_go[, -c(7)]

# Filter out GO terms related to cellular components
results_ora_go <- subset(results_ora_go, ONTOLOGY != "CC")

# Rename FDR column
colnames(results_ora_go)[7] <- "FDR"

# Filter results based on False Discovery Rate (FDR) threshold
results_ora_go <- filter(results_ora_go, as.numeric(FDR) <= 0.1)

### Remove some GO terms that refer to human diseases thus are not relevant in the context of this work
drop_GO <- !(results_ora_go$ID %in% c("GO:0043278", "GO:0072347", "GO:0009410", "GO:0006805", "GO:0018958", "GO:0016999", "GO:0060992", "GO:0017085", 
                                      "GO:0071466", "GO:0009404", "GO:0009635", "GO:0071478", "GO:0035094", "GO:0098586", "GO:0071549", "GO:0071313",
                                      "GO:0046677", "GO:0031000", "GO:0071548", "GO:0001878", "GO:0035634", "GO:0001562", "GO:0009620", "GO:0050832",
                                      "GO:0009636"))

# Subset final results by removing unwanted GO terms
results_ora_go_final <- subset(results_ora_go, drop_GO)

# Write final results to a CSV file
write.csv(as.data.frame(results_ora_go_final), 
          file = "GO_terms_enrichment_analysis/results/clusterProfiler_GO_ORA_satellite_vs_independent_Blo_.csv", 
          row.names = FALSE)
