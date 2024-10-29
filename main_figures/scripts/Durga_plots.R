# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory 
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("gridExtra")
library("grid")
library("scales")
library("ggpubr")
library("Durga")
library("svglite")

# Read the hormone data from a CSV file and store it in 'data' variable
data <- read.csv("data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv", header = TRUE)

# Convert 'morph' column to a factor with specific order
data$morph <- factor(data$morph, levels = c("Ind", "Sat", "Fae"))

# Scale down hormone values for readability
data$blo_T <- data$blo_T/1000
data$blo_A4 <- data$blo_A4/1000
data$gon_T <- data$gon_T/1000
data$gon_A4 <- data$gon_A4/1000

# Calculate hormone ratios and add them as new columns in 'data'
data <- data %>% mutate(gon_A4overT = (gon_A4/gon_T))
data <- data %>% mutate(blo_A4overT = (blo_A4/blo_T))
data <- data %>% mutate(gon_ToverA4 = (gon_T/gon_A4))
data <- data %>% mutate(blo_ToverA4 = (blo_T/blo_A4))
data <- data %>% mutate(blo_T_over_gon_T = (blo_T/gon_T))

# Filter out rows with missing values for 'blo_A4overT'
data_A4 <- data %>% filter(!is.na(blo_A4overT))

# Define custom colors for plots
my_colors <- c("#0000FF","#D02090","#FFA500")

# Calculate Cohen's d effect sizes
d11 <- DurgaDiff(data, data.col = 10, group.col = 2,
                 groups = c("Ind" = "Ind", "Sat" = "Sat", "Fae" = "Fae"),
                 effect.type = "cohens d", na.rm = TRUE)

d12 <- DurgaDiff(data, data.col = 8, group.col = 2,
                 groups = c("Ind" = "Ind", "Sat" = "Sat", "Fae" = "Fae"),
                 effect.type = "cohens d", na.rm = TRUE)

# Export the plots as SVG
svglite("main_figures/plots/Durga_plots.svg", 
    width = 8,
    height = 8, 
    pointsize = 16,
    bg = "white")

# Set up the plot layout
layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))

# Set plot parameters
par(mgp = c(4, 1, 0), mar = c(6, 6, 4, 2) + 0.05, lwd = 1.25)

# Generate Durga plot for gonadal T
p11 <- DurgaPlot(d11,
                 contrasts = "Ind - Sat, Sat - Fae, Ind - Fae",
                 box = DurgaTransparent(my_colors, 0.0),
                 box.fill = DurgaTransparent(my_colors, 0.4),
                 box.outline = FALSE,
                 box.params = list(boxwex = 0.2),
                 violin.shape = "full",
                 violin = DurgaTransparent(my_colors, 0),
                 violin.fill = DurgaTransparent(my_colors, 0.8),
                 violin.width = 0.30,
                 points = DurgaTransparent(my_colors, 0.6),
                 points.params = list(cex = 1),
                 las = 1,
                 x.axis = TRUE,
                 ef.size.violin = DurgaTransparent("blue", 0.6),
                 ef.size.violin.shape = "full",
                 ef.size = TRUE,
                 ef.size.params = list(las = 1, 
                                       cex = 1, 
                                       mgp = c(4, 1, 0)),
                 ef.size.dx = c(0, 0, 2),
                 ef.size.pch = 15,
                 ef.size.label = "Cohen's d",
                 left.ylab = "Testicular T (ng mg-1)",
                 cex.lab = 1,
                 cex.axis = 1
)

# Generate Durga plot for gonadal A4
p12 <- DurgaPlot(d12,
                 contrasts = "Ind - Sat, Sat - Fae, Ind - Fae",
                 box = DurgaTransparent(my_colors, 0.0),
                 box.fill = DurgaTransparent(my_colors, 0.4),
                 box.outline = FALSE,
                 box.params = list(boxwex = 0.2),
                 violin.shape = "full",
                 violin = DurgaTransparent(my_colors, 0),
                 violin.fill = DurgaTransparent(my_colors, 0.8),
                 violin.width = 0.30,
                 points = DurgaTransparent(my_colors, 0.6),
                 points.params = list(cex = 1),
                 las = 1,
                 x.axis = TRUE,
                 ef.size.violin = DurgaTransparent("blue", 0.6),
                 ef.size.violin.shape = "full",
                 ef.size = TRUE,
                 ef.size.params = list(las = 1, 
                                       cex = 1, 
                                       mgp = c(3.5, 1, 0)),
                 ef.size.dx = c(0, 0, 2),
                 ef.size.pch = 15,
                 ef.size.label = "Cohen's d",
                 left.ylab = "Testicular A4 (ng mg-1)",
                 cex.lab = 1,
                 cex.axis = 1,
)

# Turn off SVG plotting device
dev.off()
