# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("ggsignif")

# Read the normalized gene expression data for blood tissue
norm_counts <- read.csv(file = "DESeq2_gene_expression_analysis/normalized_counts/normalized_counts_Blood.csv", header = TRUE, row.names = 1)

# Transpose the data for easier manipulation
norm_counts_t <- t(norm_counts)
norm_counts_t <- as.data.frame(norm_counts_t)

# Add a column 'Morph' based on the sample IDs to denote the ruff morph type
norm_counts_t$Morph <- ifelse(grepl("_IND_", rownames(norm_counts_t)), "Independent",
                              ifelse(grepl("_SAT_", rownames(norm_counts_t)), "Satellite", "Faeder"))

# Rename row names to 'Sample_ID'
norm_counts_t <- tibble::rownames_to_column(norm_counts_t, "Sample_ID")

# Reshape the data into long format for better visualization
norm_counts_t_genes <- norm_counts_t %>% pivot_longer(cols = c(HSD17B2), 
                                                      names_to = "Gene",
                                                      values_to = "Expression")

# Remove unnecessary columns
norm_counts_t_genes <- norm_counts_t_genes[, -c(2:11141)]

# Convert 'Morph' column to a factor with specific order
norm_counts_t_genes$Morph <- factor(norm_counts_t_genes$Morph, levels = c("Independent", "Satellite", "Faeder"))

# Convert expression values to numeric and round them to integers
norm_counts_t_genes$Expression <- as.numeric(norm_counts_t_genes$Expression)
norm_counts_t_genes$Expression <- round(norm_counts_t_genes$Expression, digits = 0)

# Define the label for the y-axis
y_label <- expression(paste("Blood normalized gene expression"))

# Create ggplot object
p <- ggplot(data = norm_counts_t_genes,
            aes(x = Gene, 
                y = Expression)) +
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.3,
               linewidth = 0.4,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               size = 0.2,
               width = 0.8,
               linewidth = 0.4,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.4, width = 0.5, size = 0.3, alpha = 1, position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.4, width = 0.5, size = 0.3, alpha = 1, position = position_dodge(0.9)) +
  geom_point(aes(fill = Morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = norm_counts_t_genes,
              aes(y_position = 3200,
                  xmin = 0.6, 
                  xmax = 1.35,
                  annotations = "**"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.4,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = norm_counts_t_genes,
              aes(y_position = 3050, 
                  xmin = 0.6, 
                  xmax = 1,
                  annotations = "**"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.4,
              textsize = 7,
              inherit.aes = TRUE) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000), 
                     limits = c(0, 3200), 
                     labels = c("0", "1000", "2000", "3000")) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5.5))) +
  labs(y = y_label) +
  facet_grid(~ Gene, scales = "free_x") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", family = "Arial"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14.5, family = "Arial"),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = alpha("white", 1)),
        legend.title = element_blank(),
        legend.text = element_text(size = 17, family = "Arial"),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.box.just = "left",
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box.spacing = unit(0.3, "lines"),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 18, family = "Arial", face = "italic"))
        
# Save plot as SVG file
ggsave(filename = "main_figures/plots/HSD17B2_gene_expression_blood.svg", plot = p, 
       width = 2, height = 4, bg = "white", dpi = "print", device = "svg")
