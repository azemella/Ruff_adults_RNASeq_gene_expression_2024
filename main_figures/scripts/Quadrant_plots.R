# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Load required libraries
library("tidyverse")
library("ggrepel")
library("ggpubr")

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Read DESeq2 results for Faeder vs. Independent phenotype
deseq2_results <- read.csv(
  "DESeq2_gene_expression_analysis/DESeq2_results/phenotype/full_results/DGE_DESeq2_faeder_vs_independent_Phenotype_FR.csv")

# Read allele-specific expression results for Faeder phenotype
ase_results <- read.csv(
  "allele_specific_expression_analysis/other_files/Mean_ASE_11_tissues_Faeder.csv", row.names = 1)
ase_results <- ase_results %>% mutate(ASE_ratio = mean_ase * 100)
ase_results <- ase_results[, c(1, 14)]
colnames(ase_results)[1] <- "Gene_symbol"

# Read list of inversion genes
inv_genes <- read.csv("metadata/NCBI_inversion_genes_list.csv")
inv_genes <- inv_genes[, c(1)]

# Filter DESeq2 results for genes present in the inversion genes list
deseq2_results_f <- subset(deseq2_results, Gene_symbol %in% inv_genes)

# Merge DESeq2 and ASE results based on gene symbols
data <- deseq2_results_f[, c(1, 3)]
data <- merge(data, ase_results, by = "Gene_symbol")

# Define midpoint for plot annotations
x_mid = 50
y_mid = 0

# Create plot for Faeder vs. Independent phenotype
plot_fae <- data %>% 
              ggplot(aes(x = ASE_ratio, y = log2FoldChange)) +
              geom_vline(xintercept = x_mid, 
                         linewidth = 0.5, 
                         color = "black",
                         alpha = 0.5,
                         linetype = "dashed") + 
              geom_hline(yintercept = y_mid, 
                         linewidth = 0.5, 
                         color = "black",
                         alpha = 0.5,
                         linetype = "dashed") + 
              geom_text_repel(data = subset(data, Gene_symbol == "HSD17B2"),
                              aes(label = paste(expression(italic("HSD17B2"))),
                                  color = "black"), 
                              color = "black",
                              box.padding = 0.5,
                              nudge_y = 0.3,
                              nudge_x = 0,
                              size = 8,
                              parse = TRUE,
                              segment.color = 'transparent') +
              geom_text_repel(data = subset(data, Gene_symbol == "DPEP1"),
                              aes(label = paste(expression(italic("DPEP1"))),
                                  color = "black"), 
                              color = "black",
                              box.padding = 0.5,
                              nudge_y = 0.3,
                              nudge_x = 0,
                              size = 8,
                              parse = TRUE,
                              segment.color = 'transparent') +
              geom_point(aes(fill = ASE_ratio), 
                         size = 7, 
                         alpha = 1, 
                         color = "black", 
                         shape = 21) +
              scale_color_gradient2(low = "#33499B", mid = "#EFE7BB", high = "#A50026", 
                                    midpoint = 50, 
                                    limits = c(0, 100),
                                    guide = guide_colorbar(reverse = FALSE, order = 1, ticks.colour = "black", frame.colour = "black", 
                                               frame.linewidth = 1/.pt, ticks.linewidth = 1/.pt)) +
              scale_fill_gradient2(low = "#33499B", mid = "#EFE7BB", high = "#A50026", 
                                   midpoint = 50, 
                                   limits = c(0, 100),
                                   name = "ASE",
                                   guide = guide_colorbar(reverse = FALSE, order = 1, ticks.colour = "black", frame.colour = "black", 
                                              frame.linewidth = 1/.pt, ticks.linewidth = 1/.pt)) +
              scale_y_continuous(limits = c(-3, 3), 
                                 breaks = c(-3, -2, -1, 0, 1, 2, 3), 
                                 labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
              xlim(0, 100) +
              labs(x = "Faeder inversion allele proportion (%)",
                   y = expression(Biallelic~expression~(Log[2]~FC)),
                   title = "Faeder vs. Independent") +
              theme_bw() +
              theme(axis.text = element_text(size = 18, color = "black", family = "Arial"),
                    axis.title.x = element_text(size = 28, family = "Arial", margin = margin(10, 0, 0, 0)),
                    axis.title.y = element_text(size = 28, family = "Arial", margin = margin(0, 10, 0, 0)),
                    axis.ticks.length = unit(.15, "cm"),
                    axis.ticks = element_line(color = "black", linewidth = 0.5),
                    panel.border = element_rect(color = "black", fill = NA, size = 1.25),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(size = 32, hjust = 0.5, family = "Arial"), 
                    legend.position = "right",
                    legend.title = element_text(size = 22, family = "Arial"),
                    legend.text = element_text(size = 15, family = "Arial"),
                    legend.spacing.x = unit(0.2, "cm"),
                    legend.key.height = unit(1.2, "cm"),
                    aspect.ratio = 1)


# Save the plot as SVG file
ggsave(filename = "main_figures/plots/Quadrants_fae.svg", 
        plot = plot_fae, width = 8, height = 8, bg = "white", dpi = "print", device = "svg")

# Read DESeq2 results for Satellite vs. Independent phenotype
deseq2_results <- read.csv(
  "DESeq2_gene_expression_analysis/DESeq2_results/phenotype/full_results/DGE_DESeq2_satellite_vs_independent_Phenotype_FR.csv")

# Read allele-specific expression results for Satellite phenotype
ase_results <- read.csv(
  "allele_specific_expression_analysis/other_files/Mean_ASE_11_tissues_Satellite.csv", row.names = 1)
ase_results <- ase_results %>% mutate(ASE_ratio = mean_ase * 100)
ase_results <- ase_results[, c(1, 14)]
colnames(ase_results)[1] <- "Gene_symbol"

# Filter DESeq2 results for genes present in the inversion genes list
deseq2_results_f <- subset(deseq2_results, Gene_symbol %in% inv_genes)

# Merge DESeq2 and ASE results based on gene symbols
data <- deseq2_results_f[, c(1, 3)]
data <- merge(data, ase_results, by = "Gene_symbol")

# Adjust log2FoldChange values greater than 3 to 3
data$log2FoldChange <- ifelse(data$log2FoldChange > 3, 3, data$log2FoldChange)

# Define midpoint for plot annotations
x_mid = 50
y_mid = 0

# Create plot for Satellite vs. Independent phenotype
plot_sat <- data %>% 
              ggplot(aes(x = ASE_ratio, y = log2FoldChange)) +
              geom_vline(xintercept = x_mid, 
                         linewidth = 0.5, 
                         color = "black",
                         alpha = 0.5,
                         linetype = "dashed") + 
              geom_hline(yintercept = y_mid, 
                         linewidth = 0.5, 
                         color = "black",
                         alpha = 0.5,
                         linetype = "dashed") + 
              geom_text_repel(data = subset(data, Gene_symbol == "HSD17B2"),
                              aes(label = paste(expression(italic("HSD17B2"))),
                                  color = "black"), 
                              color = "black",
                              box.padding = 0.5,
                              nudge_y = 0.3,
                              nudge_x = -1.5,
                              size = 8,
                              parse = TRUE,
                              segment.color = 'transparent') +
               geom_text_repel(data = subset(data, Gene_symbol == "DPEP1"),
                               aes(label = paste(expression(italic("DPEP1"))),
                                   color = "black"), 
                               color = "black",
                               box.padding = 0.5,
                               nudge_y = -0.3,
                               nudge_x = 1.5,
                               size = 8,
                               parse = TRUE,
                               segment.color = 'transparent') +
              geom_point(data = subset(data, Gene_symbol != "TUBB3"),
                         aes(fill = ASE_ratio), 
                         size = 7, 
                         alpha = 1, 
                         color = "black", 
                         shape = 21) +
              geom_text_repel(data = subset(data, Gene_symbol == "TUBB3"),
                              aes(label = paste(expression(italic("TUBB3"))),
                                  color = "black"), 
                              color = "black",
                              box.padding = 0.5,
                              nudge_y = -0.1,
                              nudge_x = -1,
                              size = 8,
                              parse = TRUE,
                              segment.color = 'transparent') +
              geom_point(data = subset(data, Gene_symbol == "TUBB3"),
                         aes(fill = ASE_ratio), 
                         size = 7, 
                         alpha = 1, 
                         color = "black", 
                         shape = 24) +
              scale_color_gradient2(low = "#33499B", mid = "#EFE7BB", high = "#A50026", 
                                    midpoint = 50, 
                                    limits = c(0, 100),
                                    guide = guide_colorbar(reverse = FALSE, order = 1, ticks.colour = "black", frame.colour = "black", 
                                               frame.linewidth = 1/.pt, ticks.linewidth = 1/.pt)) +
              scale_fill_gradient2(low = "#33499B", mid = "#EFE7BB", high = "#A50026", 
                                   midpoint = 50, 
                                   limits = c(0, 100),
                                   name = "ASE",
                                   guide = guide_colorbar(reverse = FALSE, order = 1, ticks.colour = "black", frame.colour = "black", 
                                              frame.linewidth = 1/.pt, ticks.linewidth = 1/.pt)) +
              scale_y_continuous(limits = c(-3, 3), 
                                 breaks = c(-3, -2, -1, 0, 1, 2, 3), 
                                 labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
              xlim(0, 100) +
              labs(x = "Satellite inversion allele proportion (%)",
                   y = expression(Biallelic~expression~(Log[2]~FC)),
                   title = "Satellite vs. Independent") +
              theme_bw() +
              theme(axis.text = element_text(size = 18, color = "black", family = "Arial"),
                    axis.title.x = element_text(size = 28, family = "Arial", margin = margin(10, 0, 0, 0)),
                    axis.title.y = element_text(size = 28, family = "Arial", margin = margin(0, 10, 0, 0)),
                    axis.ticks.length = unit(.15, "cm"),
                    axis.ticks = element_line(color = "black", linewidth = 0.5),
                    panel.border = element_rect(color = "black", fill = NA, size = 1.25),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(size = 32, hjust = 0.5, family = "Arial"), 
                    legend.position = "right",
                    legend.title = element_text(size = 22, family = "Arial"),
                    legend.text = element_text(size = 15, family = "Arial"),
                    legend.spacing.x = unit(0.2, "cm"),
                    legend.key.height = unit(1.2, "cm"),
                    aspect.ratio = 1) 

# Save the plot as SVG file
ggsave(filename = "main_figures/plots/Quadrants_sat.svg", 
        plot = plot_sat, width = 8, height = 8, bg = "white", dpi = "print", device = "svg")