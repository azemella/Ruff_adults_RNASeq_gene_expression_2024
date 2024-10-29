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

# Remove NAs
data <- data[!is.na(data$blo_A4),]

# Create ggplot object
p <- ggplot(data, aes(x = gon_A4, y = blo_A4)) +
  geom_point(aes(fill = morph, color = morph), size = 3.5) +
  geom_smooth(aes(fill = morph, color = morph), method = "lm", alpha = 0.2) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE), expand = c(0, 0), limits = c(-2.5, 23)) +
  scale_x_continuous(breaks = c(0.000, 0.025, 0.050, 0.0750, 0.100),
                     labels = c("0", "0.025", "0.05", "0.075", "0.1")) +
  labs(x = expression("Testicular A4 (ng mg"^{-1}*")"),
       y = expression("Blood A4 (ng ml"^{-1}*")")) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 10, fill = NA))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 24, colour = "black", family = "Arial"),
    axis.text.y = element_text(size = 24, colour = "black", family = "Arial"),
    axis.ticks.y = element_line(color = "black", linewidth = 0.75),
    axis.ticks.x = element_line(color = "black", linewidth = 0.75),
    axis.ticks.length = unit(.25, "cm"),
    axis.title.x.bottom = element_text(size = 36, family = "Arial",margin = margin(5, 0, 0, 0)),
    axis.title.y.left = element_text(size = 36, family = "Arial",margin = margin(0, 15, 0, 0)),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    legend.text = element_text(size = 36, family = "Arial"),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.key = element_blank()
  ) 

# Save plot as SVG file
ggsave(filename = "supplementary_figures/plots/Scatterplot_BloodA4_TestA4.svg", plot = p, 
       width = 10, height = 10, bg = "white", dpi = "print", device = "svg")
