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
data <- data[!is.na(data$blo_A4overT),]

# Define y-axis label
y_label <- expression(paste("Blood A4/T ratio"))

# Create ggplot object
p <- ggplot(data = data,
            aes(x = morph, 
                y = blo_A4overT)) +
  geom_point(size = -1, aes(fill = morph, color = morph)) +
  geom_boxplot(aes(fill = morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  labs(y = y_label) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150),
                     labels = c("0", "25", "50", "75", "100", "125", "150"),
                     limits = c(0, 150)) +
  theme(
    axis.text.x = element_text(size = 36, colour = "black", family = "Arial"),
    axis.text.y = element_text(size = 24, colour = "black", family = "Arial"),
    axis.ticks.y = element_line(color = "black", linewidth = 0.75),
    axis.ticks.x = element_line(color = "black", linewidth = 0.75),
    axis.ticks.length = unit(.25, "cm"),
    axis.title.x = element_blank(),
    axis.title.y.left = element_text(size = 36, family = "Arial",margin = margin(0, 15, 0, 0)),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white", 1)),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "none"
  ) 

# Save plot as SVG file
ggsave(filename = "supplementary_figures/plots/Boxplot_BloodA4_overT.svg", plot = p, 
       width = 10, height = 10, bg = "white", dpi = "print", device = "svg")
