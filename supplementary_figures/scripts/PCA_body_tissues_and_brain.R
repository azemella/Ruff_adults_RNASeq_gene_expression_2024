# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load necessary libraries
library("DESeq2")
library("ashr")
library("BiocParallel")
library("tidyverse")
library("limma")
library("scales")

# Read metadata and preprocess it
samples <- read.csv(file.path(dir, "metadata/metadata.csv")) %>% distinct(Sample_ID, .keep_all = TRUE)
samples <- samples[samples$Tissue != "Blo",]
samples <- samples %>% mutate(Tissue = ifelse(Tissue %in% c("A+TnA", "VTA+SN", "POM", "Hyp", "NC", "Rap", "LS"), "Bra", Tissue))
samples$Phenotype <- factor(samples$Phenotype)
samples$Tissue <- factor(samples$Tissue)
rownames(samples) <- samples$Sample_ID

# Read count matrix
count_matrix <- read.csv("DESeq2_gene_expression_analysis/STAR_gene_counts.csv", row.names = 1)
count_matrix <- count_matrix[, rownames(samples)]

# Design matrix
design <- ~ Tissue + Phenotype
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = design)

# Filter low count genes
ind_sample_count <- sum(samples$phenotype == "Independent")
sat_sample_count <- sum(samples$phenotype == "Satellite")
fae_sample_count <- sum(samples$phenotype == "Faeder")
smallestGroupSize <- min(ind_sample_count, sat_sample_count, fae_sample_count)
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

# Perform DE analysis
dds <- DESeq(dds, parallel = TRUE)
dds <- dds[which(mcols(dds)$betaConv),]

# Perform variance stabilizing transformation
vst <- varianceStabilizingTransformation(dds, blind = FALSE)

# Perform PCA analysis
pcaData <- plotPCA(vst, intgroup = c("Tissue", "Phenotype"), returnData = TRUE)
pcaData$Tissue <- factor(pcaData$Tissue, levels = c("Pit", "Gon", "Adr", "Liv", "Bra"))
pcaData$Phenotype <- factor(pcaData$Phenotype, levels = c("Independent", "Satellite", "Faeder"))
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create ggplot object
p <- ggplot(pcaData, aes(PC1, PC2, fill = Tissue, shape = Phenotype)) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    geom_point(aes(fill = Tissue), 
               size = 5.5,
               color = "black",
               stroke = .35,
               shape = 21,
               position = "jitter") +
    scale_fill_manual(name = "Tissue", 
                      labels = c("PIT", "TES", "ADR", "LIV", "BRA"),
                      values = alpha(c('midnightblue', "springgreen3", 'yellow', 'violet', 'grey'), 1)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_x_continuous(labels = label_number(accuracy = 1)) +
    theme(
      axis.text.x = element_text(size = 16, color = alpha("black", 1), family = "Arial"), 
      axis.title.x = element_text(size = 28, family = "Arial", vjust = -0.25),
      axis.text.y = element_text(size = 16, color = alpha("black", 1), family = "Arial"), 
      axis.title.y = element_text(size = 28, family = "Arial", vjust = 1.5),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks.x.bottom = element_line(linewidth = 1, color = "black"),
      axis.ticks.y.left = element_line(linewidth = 1, color = "black"),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.25),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 0.75
  ) +
    guides(fill = guide_legend(override.aes = list(shape = 21), order = 1),
           shape = guide_legend(order = 2))

# Save plot as SVG file
ggsave(filename = "supplementary_figures/plots/PCA_Brain_and_Body_Tissues.svg", plot = p, 
       width = 10, height = 8, bg = "white", dpi = "print", device = "svg")
