# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("ggplot2")
library("dplyr")
library("forcats")
library("ggsci")
library("ggsignif")
library("scales")
library("gridExtra")
library("grid")

# Read normalized counts data
norm_counts <- read.csv(file = "DESeq2_gene_expression_analysis/normalized_counts/normalized_counts_Phenotype.csv", header = TRUE, row.names = 1)
norm_counts <- norm_counts[, -c(29:37)]

# Transpose the data frame
norm_counts_t <- t(norm_counts)

# Add 'Morph' column based on sample names
if (is.data.frame(norm_counts_t)) {
  norm_counts_t$Morph <- ifelse(grepl("_I_", rownames(norm_counts_t)), "Independent",
                                ifelse(grepl("_S_", rownames(norm_counts_t)), "Satellite", "Faeder"))
} else if (is.matrix(norm_counts_t)) {
  norm_counts_t <- cbind(norm_counts_t, Morph = ifelse(grepl("_I_", rownames(norm_counts_t)), "Independent",
                                                       ifelse(grepl("_S_", rownames(norm_counts_t)), "Satellite", "Faeder")))
}

# Convert to data frame and add 'Tissue' column based on sample names
norm_counts_t <- as.data.frame(norm_counts_t)
norm_counts_t <- norm_counts_t %>%
  mutate(Tissue = case_when(
    grepl("_ADR_", rownames(norm_counts_t)) ~ "ADR",
    grepl("_GON_", rownames(norm_counts_t)) ~ "TES",
    grepl("_HYP_", rownames(norm_counts_t)) ~ "HYP",
    grepl("_LS_", rownames(norm_counts_t)) ~ "LS",
    grepl("_NC_", rownames(norm_counts_t)) ~ "NC",
    grepl("_PIT_", rownames(norm_counts_t)) ~ "PIT",
    grepl("_POA_", rownames(norm_counts_t)) ~ "POM",
    grepl("_RAP_", rownames(norm_counts_t)) ~ "RAP",
    grepl("_TNA_", rownames(norm_counts_t)) ~ "A+TnA",
    grepl("_VTA_", rownames(norm_counts_t)) ~ "VTA+SN",
    grepl("_LIV_", rownames(norm_counts_t)) ~ "LIV",
    TRUE ~ NA_character_
  ))

# Rename rownames column to "Sample_ID"
norm_counts_t <- tibble::rownames_to_column(norm_counts_t, "Sample_ID")

# Select relevant columns and rename
norm_counts_t <- norm_counts_t[, c("Sample_ID", "PGR", "Morph", "Tissue")]
norm_counts_t$Gene <- "PGR"

# Factorize 'Morph' and 'Tissue' columns
norm_counts_t$Morph <- factor(norm_counts_t$Morph, levels = c("Independent", "Satellite", "Faeder"))
norm_counts_t$Tissue <- factor(norm_counts_t$Tissue, levels = c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV"))

# Rename column to "Expression" and ensure numeric type
colnames(norm_counts_t)[2] <- "Expression"
norm_counts_t$Expression <- as.numeric(norm_counts_t$Expression)

# Adjust negative expression values to zero
norm_counts_t$Expression <- ifelse(norm_counts_t$Expression < 0, 0, norm_counts_t$Expression)

# Define y-axis label
y_label <- expression(paste(italic("PGR") ~ "normalized gene expression"))

# Create a data frame for plotting
dat <- data.frame(Tissue = c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV"))
dat$Tissue <- factor(dat$Tissue, levels = c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV"))

# Create ggplot object
p <- ggplot(data = norm_counts_t,
            aes(x = Tissue, 
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
               show.legend = FALSE) +
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
              aes(y_position = c(NA, NA, NA, NA, NA, NA, NA, 3950, NA, NA, NA), 
                  xmin = c(NA, NA, NA, NA, NA, NA, NA, 0.65, NA, NA, NA), 
                  xmax = c(NA, NA, NA, NA, NA, NA, NA, 1.35, NA, NA, NA),
                  annotations = c("", "", "", "", "", "", "", "***", "", "", "")), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 6,
              inherit.aes = TRUE) +
  geom_signif(data = dat,
              aes(y_position = c(NA, NA, NA, NA, NA, NA, NA, 3750, NA, NA, NA), 
                  xmin = c(NA, NA, NA, NA, NA, NA, NA, 1, NA, NA, NA), 
                  xmax = c(NA, NA, NA, NA, NA, NA, NA, 1.35, NA, NA, NA),
                  annotations = c("", "", "", "", "", "", "", "**", "", "", "")), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 6,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c("0", "1000", "2000", "3000", "4000"),
                     limits = c(0, 4000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  labs(y = y_label) +
  facet_grid(~ Tissue, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, family = "Arial"),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white")),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, family = "Arial"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 18, family = "Arial"),
  ) 

# Save plot as SVG file
ggsave(filename = "supplementary_figures/plots/PGR_gene_expression_all_tissues_Norm.svg", plot = p, 
       width = 10, height = 4, bg = "white", dpi = "print", device = "svg")
