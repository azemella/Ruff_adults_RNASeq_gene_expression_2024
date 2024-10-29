# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Load necessary libraries
library("DESeq2")
library("ashr")
library("tibble")
library("dplyr")
library("BiocParallel")

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Register parallel processing
register(MulticoreParam(4))

# To analyze a different tissue, subset the metadata file for a different tissue
# List of available tissues: (Adr, A+TnA, Blo, Gon, Hyp, Liv, LS, NC, Pit, POM, Rap, VTA+SN)
# Batch effect (flow cell code) has to be accounted when analyzing certain tissues
# Change the script to remove/add "Batch" from the exp design (mask or unmask lines #38-#39 accordingly)

# Read metadata, keep only distinct samples
samples <- read.csv(file.path(dir, "metadata/metadata.csv")) %>% distinct(Sample_ID, .keep_all = TRUE)

# Here we subset for Adrenal glands
samples <- samples[samples$Tissue == "Adr",]
samples$Phenotype <- factor(samples$Phenotype)
colnames(samples)[22] <- "Batch"
samples$Batch <- factor(samples$Batch)
rownames(samples) <- samples$Sample_ID

# Read count matrix
count_matrix <- read.csv("DESeq2_gene_expression_analysis/STAR_gene_counts.csv", row.names = 1)
count_matrix <- count_matrix[, rownames(samples)]

# Design matrix
design <- ~ Batch + Phenotype
#design <- ~ Phenotype
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = design)

# Filter low count genes
ind_sample_count <- sum(samples$Phenotype == "Independent")
sat_sample_count <- sum(samples$Phenotype == "Satellite")
fae_sample_count <- sum(samples$Phenotype == "Faeder")
smallestGroupSize <- min(ind_sample_count, sat_sample_count, fae_sample_count)
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

# Perform DE analysis
dds <- DESeq(dds, parallel = TRUE)
dds <- dds[which(mcols(dds)$betaConv),]

# Contrast comparisons
res1 <- results(dds, contrast = c("Phenotype", "Faeder", "Independent"), alpha = 0.05, parallel = TRUE)
res2 <- results(dds, contrast = c("Phenotype", "Satellite", "Independent"), alpha = 0.05, parallel = TRUE)
res3 <- results(dds, contrast = c("Phenotype", "Faeder", "Satellite"), alpha = 0.05, parallel = TRUE)

# Summarize results
summary(res1)
summary(res2)
summary(res3)

# Shrinkage estimation
resShrink1 <- lfcShrink(dds = dds, res = res1, type = "ashr")
resShrink2 <- lfcShrink(dds = dds, res = res2, type = "ashr")
resShrink3 <- lfcShrink(dds = dds, res = res3, type = "ashr")

# Convert results to data frames
full_res1 <- as.data.frame(resShrink1)
full_res1 <- tibble::rownames_to_column(full_res1, "Gene_symbol")

full_res2 <- as.data.frame(resShrink2)
full_res2 <- tibble::rownames_to_column(full_res2, "Gene_symbol")

full_res3 <- as.data.frame(resShrink3)
full_res3 <- tibble::rownames_to_column(full_res3, "Gene_symbol")

# Write full results to CSV files
write.csv(as.data.frame(full_res1), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/full_results/DGE_DESeq2_faeder_vs_independent_Adr_FR.csv", row.names = F)
write.csv(as.data.frame(full_res2), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/full_results/DGE_DESeq2_satellite_vs_independent_Adr_FR.csv", row.names = F)
write.csv(as.data.frame(full_res3), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/full_results/DGE_DESeq2_faeder_vs_satellite_Adr_FR.csv", row.names = F)

# Subset significant results
sig_res1 <- as.data.frame(resShrink1)
sig_res1 <- tibble::rownames_to_column(sig_res1, "Gene_symbol")
sig_res1 <- sig_res1[order(sig_res1$padj),]
sig_res1 <- subset(sig_res1, padj < 0.05)

sig_res2 <- as.data.frame(resShrink2)
sig_res2 <- tibble::rownames_to_column(sig_res2, "Gene_symbol")
sig_res2 <- sig_res2[order(sig_res2$padj),]
sig_res2 <- subset(sig_res2, padj < 0.05)

sig_res3 <- as.data.frame(resShrink3)
sig_res3 <- tibble::rownames_to_column(sig_res3, "Gene_symbol")
sig_res3 <- sig_res3[order(sig_res3$padj),]
sig_res3 <- subset(sig_res3, padj < 0.05)

# Write significant results to CSV files
write.csv(as.data.frame(sig_res1), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/significant_results/DGE_DESeq2_faeder_vs_independent_Adr_SR.csv", row.names = F)
write.csv(as.data.frame(sig_res2), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/significant_results/DGE_DESeq2_satellite_vs_independent_Adr_SR.csv", row.names = F)
write.csv(as.data.frame(sig_res3), 
          file = "DESeq2_gene_expression_analysis/DESeq2_results/tissue/significant_results/DGE_DESeq2_faeder_vs_satellite_Adr_SR.csv", row.names = F)
