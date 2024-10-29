# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Loading required libraries
library("tidyverse")

# Read the data files containing allele-specific expression (ASE) counts for different morphs
ase_counts_fae <- read.csv(file = "others/HSD17B2_ASE_Faeder.csv", header = TRUE)
ase_counts_sat <- read.csv(file = "others/HSD17B2_ASE_Satellite.csv", header = TRUE)

# Combine the ASE counts from both inversion morphs into a single dataset
ase_counts <- rbind(ase_counts_fae, ase_counts_sat)

# Filter the data to include only samples from blood tissue
ase_counts <- ase_counts[ase_counts$Tissue == "Blood", ]

# Rename column names
colnames(ase_counts)[1] <- "sample_name"
colnames(ase_counts)[2] <- "inversion_allele_ratio"

# Add a new column to denote the gene being analyzed (HSD17B2)
ase_counts$Gene <- "HSD17B2"

# Convert the 'Morph' column to a factor and set the order
ase_counts$Morph <- factor(ase_counts$Morph, levels = c("Satellite", "Faeder"))

# Define the label for the y-axis
y_label <- expression("Blood inversion allele proportion (%)")

# Convert inversion allele ratios to percentages and round up
ase_counts$inversion_allele_ratio <- ase_counts$inversion_allele_ratio * 100
ase_counts$inversion_allele_ratio <- ceiling(ase_counts$inversion_allele_ratio)

# Create the ggplot object
p <- ggplot(data = ase_counts,
            aes(x = as.factor(Tissue), 
                y = inversion_allele_ratio)) +
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_hline(yintercept = 50, 
             color = "black",
             alpha = 0.5,
             linewidth = 0.75,
             linetype = "dashed",
             show.legend = FALSE) +
  geom_boxplot(aes(fill = Morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.4,
               outlier.shape = 21,
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
  scale_fill_manual(labels = c("Satellite","Faeder"), 
                    values = c("violetred", "orange")) +
  scale_color_manual(labels = c("Satellite","Faeder"), 
                     values = c("violetred", "orange")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 5, linetype = 0))) +
  labs(y = y_label) +
  ylim(0, 100) +
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
ggsave(filename = "main_figures/plots/HSD17B2_allelic_imbalance_blood.svg", plot = p, 
       width = 2, height = 4, bg = "white", dpi = "print", device = "svg")
