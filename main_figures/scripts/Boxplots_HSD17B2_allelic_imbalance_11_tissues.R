# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("plotly")
library("svglite")

# Read the allele-specific expression (ASE) counts for Faeder and Satellite morphs
ase_counts_fae <- read.csv(file = "others/HSD17B2_ASE_Faeder.csv", header = TRUE)
ase_counts_sat <- read.csv(file = "others/HSD17B2_ASE_Satellite.csv", header = TRUE)

# Combine the ASE counts from Faeder and Satellite morphs into a single dataset
ase_counts <- rbind(ase_counts_fae, ase_counts_sat)

# Rename column names
colnames(ase_counts)[1] <- "sample_name"
colnames(ase_counts)[2] <- "inversion_allele_ratio"
ase_counts$Gene <- "HSD17B2"

# Update tissue abbreviations
ase_counts$Tissue[ase_counts$Tissue == "GON"] <- "TES"
ase_counts$Tissue[ase_counts$Tissue == "POA"] <- "POM"

# Convert 'Morph' and 'Tissue' columns to factors with specific order
ase_counts$Morph <- factor(ase_counts$Morph, levels = c("Satellite", "Faeder"))
ase_counts$Tissue <- factor(ase_counts$Tissue, levels = c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV"))

# Define the label for the y-axis
y_label <- expression(paste(italic("HSD17B2"), " inversion allele proportion (%)"))

# Convert inversion allele ratios to percentages and round them up
ase_counts$inversion_allele_ratio <- ase_counts$inversion_allele_ratio * 100
ase_counts$inversion_allele_ratio <- ceiling(ase_counts$inversion_allele_ratio)

# Remove rows with missing values
ase_counts <- ase_counts[complete.cases(ase_counts), ]

# Create ggplot object
p <- ggplot(data = ase_counts,
            aes(x = as.factor(Tissue), 
                y = inversion_allele_ratio)) +
        geom_point(size = -1, aes(fill = Morph, color = Morph)) +
        geom_hline(yintercept = 50, 
                   color = "black",
                   alpha = 0.5,
                   linewidth = 0.4,
                   linetype = "dashed",
                   show.legend = FALSE) +
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
        scale_fill_manual(labels = c("Satellite","Faeder"), 
                          values = c("violetred", "orange")) +
        scale_color_manual(labels = c("Satellite","Faeder"), 
                           values = c("violetred", "orange")) +
        guides(color = guide_legend(override.aes = list(shape = 21, size = 4, linetype = 0))) +
        labs(y = y_label) +
        ylim(0, 100) +
        facet_grid(~ Tissue, scales = "free_x") +
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
              strip.text.x = element_text(size = 18, family = "Arial"))

# Save plot as SVG file
ggsave(filename = "main_figures/plots/HSD17B2_allelic_imbalance_tissues.svg", plot = p, 
       width = 10, height = 4, bg = "white", dpi = "print", device = "svg")
