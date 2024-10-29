# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("ggrepel")

# Read DESeq2 results for Faeder vs. Satellite phenotypes
deseq2_results <- read.csv("DESeq2_gene_expression_analysis/DESeq2_results/tissue/full_results/DGE_DESeq2_faeder_vs_satellite_Hyp_FR.csv")

# Calculate -log10 adjusted p-values
deseq2_results$log10padj <- -log10(deseq2_results$padj)

# Cap the maximum value of log10padj at 20
deseq2_results$log10padj <- ifelse(deseq2_results$log10padj > 20, 20, deseq2_results$log10padj)

# Cap log2FoldChange values at -3 and 3
deseq2_results$log2FoldChange <- ifelse(deseq2_results$log2FoldChange > 3, 3,
                                        ifelse(deseq2_results$log2FoldChange < -3, -3, deseq2_results$log2FoldChange))

# Assign colors based on significance and fold change
colors <- ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0.2, "Upregulated",
                 ifelse(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < -0.2, "Downregulated", "NS"))

# Assign shapes based on significance and fold change
shapes <- ifelse(deseq2_results$log10padj == 20 | deseq2_results$log2FoldChange == 3 | deseq2_results$log2FoldChange == -3, "Triangle", "Circle")

# Select top differentially expressed inversion genes for highlighting
genes <- deseq2_results[deseq2_results$Gene_symbol %in% c("PLCG2", "GINS2", "TUBB3", "HSD17B2", "SLC7A5"), ]

# Create ggplot object
p <- ggplot(deseq2_results, aes(x = log2FoldChange, y = log10padj, color = colors, shape = shapes, fontface = "italic")) + 
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% c("HSD17B2", "SLC7A5", "PLCG2", "GINS2"), ],
             aes(x = log2FoldChange, y = log10padj), 
             size = 6, color = "black", shape = 19) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% "TUBB3", ],
             aes(x = log2FoldChange, y = log10padj), 
             size = 6, color = "black", shape = 17) +
  geom_point(size = 5) +
  scale_color_manual(values = alpha(c("Upregulated" = "red", "Downregulated" = "blue", "NS" = "darkgrey"), .35),
                     breaks = c("Downregulated", "Upregulated", "NS")) +
  scale_shape_manual(values = c("Triangle" = 17, "Circle" = 19), guide = "none") +
  xlab(bquote("Log"[2] ~ "fold change")) +
  ylab(bquote("-Log"[10] ~ italic("P") * "-adjust")) +
  labs(color = "Gene significance") +
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 10, 20)) +
  geom_vline(xintercept = c(0.2, -0.2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.25, linetype = "dashed", color = "black") +
  ggtitle("Faeder vs. Satellite") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 22, color = alpha("black", 1)),
    axis.text.y = element_text(size = 22, color = alpha("black", 1)),
    axis.title.x = element_text(size = 36, color = alpha("black", 1), vjust = -0.15),
    axis.title.y = element_text(size = 36, color = alpha("black", 1), vjust = 1),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks.x.bottom = element_line(linewidth = 1, color = "black"),
    axis.ticks.y.left = element_line(linewidth = 1, color = "black"),
    legend.position = "none",
    plot.title = element_text(size = 36, color = alpha("black", 1), hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  geom_text_repel(data = deseq2_results[deseq2_results$Gene_symbol %in% "SLC7A5", ],
                  aes(x = log2FoldChange, 
                      y = log10padj, 
                      label = "SLC7A5",
                      color = "black", 
                      shape = "Triangle"), 
                  color = "black",
                  box.padding = 0.5,
                  nudge_y = 1,
                  nudge_x = 0.75,
                  size = 8) +
  geom_text_repel(data = deseq2_results[deseq2_results$Gene_symbol %in% "HSD17B2", ],
                  aes(x = log2FoldChange, 
                      y = log10padj, 
                      label = "HSD17B2",
                      color = "black", 
                      shape = "Triangle"), 
                  color = "black",
                  box.padding = 0.5,
                  nudge_y = 0.75,
                  nudge_x = 0.75,
                  size = 8) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% c("HSD17B2", "SLC7A5"), ],
             aes(x = log2FoldChange, y = log10padj), 
             size = 6, color = "black", shape = 19) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% c("HSD17B2", "SLC7A5"), ],
             size = 5, color = "red", shape = 19) +
  geom_text_repel(data = deseq2_results[deseq2_results$Gene_symbol %in% "TUBB3", ],
                  aes(x = log2FoldChange, 
                      y = log10padj, 
                      label = "TUBB3",
                      color = "black", 
                      shape = "Triangle"), 
                  color = "black",
                  box.padding = 0.5,
                  nudge_y = -0.5,
                  nudge_x = 0.75,
                  size = 8) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% "TUBB3", ],
             aes(x = log2FoldChange, y = log10padj), 
             size = 6, color = "black", shape = 17) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% "TUBB3", ],
             size = 5, color = "blue", shape = 17) +
  geom_text_repel(data = deseq2_results[deseq2_results$Gene_symbol %in% "GINS2", ],
                  aes(x = log2FoldChange, 
                      y = log10padj, 
                      label = "GINS2",
                      color = "black", 
                      shape = "Triangle"), 
                  color = "black",
                  box.padding = 0.5,
                  nudge_y = 0.75,
                  nudge_x = -0.75,
                  size = 8) +
  geom_text_repel(data = deseq2_results[deseq2_results$Gene_symbol %in% "PLCG2", ],
                  aes(x = log2FoldChange, 
                      y = log10padj, 
                      label = "PLCG2",
                      color = "black", 
                      shape = "Triangle"), 
                  color = "black",
                  box.padding = 0.5,
                  nudge_y = 0.75,
                  nudge_x = -0.75,
                  size = 8) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% c("GINS2", "PLCG2"), ],
             aes(x = log2FoldChange, y = log10padj), 
             size = 6, color = "black", shape = 19) +
  geom_point(data = deseq2_results[deseq2_results$Gene_symbol %in% c("GINS2", "PLCG2"), ],
             size = 5, color = "blue", shape = 19)

# Save plot as SVG file
ggsave(filename = "main_figures/plots/VolcanoPlot_Hyp_Fae_vs_Sat.svg", plot = p, 
       width = 8, height = 8, bg = "white", dpi = "print", device = "svg")
