# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tibble")
library("dplyr")

# Read the list of differentially expressed genes
sign_genes <- read.csv("DESeq2_gene_expression_analysis/DESeq2_results/phenotype/significant_results/DGE_DESeq2_faeder_vs_satellite_Phenotype_SR.csv")
sign_genes <- sign_genes[, 1]

# Read the list of all background genes
all_genes <- read.csv("DESeq2_gene_expression_analysis/DESeq2_results/phenotype/full_results/DGE_DESeq2_faeder_vs_satellite_Phenotype_FR.csv")
all_genes <- all_genes[, 1]

# Read the list of genes within the inversion region boundaries and filter out those present in the list of all genes
inversion_genes <- read.csv("metadata/NCBI_inversion_genes_list.csv")
inversion_genes <- inversion_genes[, 2]
inversion_genes <- intersect(inversion_genes, all_genes)

# Calculate the number of genes that are both differentially expressed and within the inversion region
a <- intersect(sign_genes, inversion_genes)
a <- length(a)

# Calculate the number of genes that are differentially expressed but not within the inversion region
b <- length(sign_genes) - a

# Calculate the number of genes within the inversion region but not differentially expressed
c <- length(intersect(inversion_genes, all_genes)) - a

# Calculate the number of genes that are neither differentially expressed nor within the inversion region
d <- length(all_genes) - length(sign_genes) - length(inversion_genes) - length(intersect(inversion_genes, sign_genes))

# Create a 2x2 contingency table
dat <- data.frame(
  "Differentially_expressed" = c(a, b),
  "Not_differentially_expressed" = c(c, d),
  row.names = c("Inside_inversion", "Outside_inversion"),
  stringsAsFactors = FALSE
)

# Display the data frame
dat

# Perform Fisher's exact test
test <- fisher.test(dat)
test

# Save results for all morph-tissue comparisons in a data frame and name it "Stats_enrichment_inversion_genes_DGE.csv"
# Then, run the following lines to correct p-values for multiple comparisons

# Read the table containing statistical analyses results for enrichment in inversion genes
table <- read.csv("statistical_analyses/Stats_enrichment_inversion_genes_DGE.csv")

# Adjust the p-values using the Benjamini-Hochberg correction method
table$Padjust <- p.adjust(table$Pvalue, method = "BH")

# Move the adjusted p-values next to the original p-values in the table
table <- table %>% relocate(Padjust, .after = Pvalue)

# Save results
write.csv(as.data.frame(table), 
          file = "statistical_analyses/Stats_enrichment_inversion_genes_DGE.csv", row.names = FALSE)

