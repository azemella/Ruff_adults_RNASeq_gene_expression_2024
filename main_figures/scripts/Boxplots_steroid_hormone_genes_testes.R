# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("ggsignif")
library("scales")
library("ggh4x")

# Read the normalized gene expression data for gonadal tissue
norm_counts <- read.csv(file = "DESeq2_gene_expression_analysis/normalized_counts/normalized_counts_Gon.csv", header = TRUE, row.names = 1)

# Select columns containing '_GON_' in their names
norm_counts <- norm_counts %>% select(contains("_GON_"))

# Transpose the data
norm_counts_t <- t(norm_counts)
norm_counts_t <- as.data.frame(norm_counts_t)

# Add a column 'Morph' based on the sample IDs 
norm_counts_t$Morph <- ifelse(grepl("_I_", rownames(norm_counts_t)), "Independent",
                              ifelse(grepl("_S_", rownames(norm_counts_t)), "Satellite", "Faeder"))

# Rename row names to 'Sample_ID'
norm_counts_t <- tibble::rownames_to_column(norm_counts_t, "Sample_ID")

# Reshape the data into long format
norm_counts_t_genes <- norm_counts_t %>% pivot_longer(cols = c(STAR, LOC106897604, LOC106896584, HSD17B3, LOC106896857), 
                                                      names_to = "Gene",
                                                      values_to = "Expression")

# Remove unnecessary columns
norm_counts_t_genes <- norm_counts_t_genes[, -c(2:18576)]

# Convert 'Morph' column to a factor with specific order
norm_counts_t_genes$Morph <- factor(norm_counts_t_genes$Morph, levels = c("Independent", "Satellite", "Faeder"))

# Convert expression values to numeric and round them to integers
norm_counts_t_genes$Expression <- as.numeric(norm_counts_t_genes$Expression)
norm_counts_t_genes$Expression <- round(norm_counts_t_genes$Expression, digits = 0)

# Update gene names
norm_counts_t_genes$Gene[norm_counts_t_genes$Gene == "LOC106897604"] <- "CYP17A1"
norm_counts_t_genes$Gene[norm_counts_t_genes$Gene == "LOC106896584"] <- "CYP11A1"
norm_counts_t_genes$Gene[norm_counts_t_genes$Gene == "LOC106896857"] <- "HSD3B2"
norm_counts_t_genes$Gene <- factor(norm_counts_t_genes$Gene, levels = c("STAR", "CYP11A1", "CYP17A1", "HSD3B2", "HSD17B3"))

# Create a data frame for gene order in plots
dat <- data.frame(Gene = c("STAR", "CYP11A1", "CYP17A1", "HSD3B2", "HSD17B3"))
dat$Gene <- factor(dat$Gene, levels = c("STAR", "CYP11A1", "CYP17A1", "HSD3B2", "HSD17B3"))

# Define the label for the y-axis
y_label <- expression(paste("Testicular normalized gene expression"))

# Create ggplot object
p <- ggplot(data = norm_counts_t_genes,
            aes(x = Gene, 
                y = Expression)) +
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.4,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               width = 0.8,
               linewidth = 0.4,
               position = position_dodge(0.9),
               show.legend = FALSE,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.4, width = 0.5,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.4, width = 0.5,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = dat,
              aes(y_position = c(7900, 2955, 39500, NA, NA), 
                  xmin = c(0.65, 0.65, 0.65, NA, NA), 
                  xmax = c(1.35, 1.35, 1.35, NA, NA),
                  annotations = c("***", "***", "***", "", "")), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = dat,
              aes(y_position = c(7450, NA, 37350, NA, NA), 
                  xmin = c(0.65, NA, 0.65, NA, NA), 
                  xmax = c(1, NA, 1, NA, NA),
                  annotations = c("***", "", "*", "", "")), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5.5))) +
  labs(y = y_label) +
  ggh4x::facet_grid2(.~Gene, scales = "free", independent = "y") +
  ggh4x::facetted_pos_scales(
    y = list(
      Gene == "STAR" ~ scale_y_continuous(breaks = c(0, 2000, 4000, 6000, 8000), limits = c(0, 8000)),
      Gene == "CYP11A1" ~ scale_y_continuous(breaks = c(0, 1000, 2000, 3000), limits = c(0, 3000)),
      Gene == "CYP17A1" ~ scale_y_continuous(breaks = c(0, 10000, 20000, 30000, 40000), limits = c(0, 40000)),
      Gene == "HSD3B2" ~ scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000), limits = c(0, 20000)),
      Gene == "HSD17B3" ~ scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000), limits = c(0, 4000)))) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14.5, family = "Arial"),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white")),
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
    strip.text.x = element_text(size = 18, family = "Arial", face = "italic"),
  ) 

# Save plot as SVG file
ggsave(filename = "Boxplots_steroid_hormone_genes_testes.svg", plot = p, 
       width = 8, height = 4, bg = "white", dpi = "print", device = "svg")
