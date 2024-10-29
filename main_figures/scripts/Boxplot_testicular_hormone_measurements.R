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
library("ggsignif")

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

# Create ggplot object
p1 <- data %>%
  ggplot(aes(fill = morph, y = gon_T, x = morph)) + 
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
             size = 5, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  geom_signif(data = data,
              aes(y_position = 0.475, 
                  xmin = 0.85, 
                  xmax = 2.15,
                  annotations = "**"), 
              tip_length = 0.02, 
              manual = T,
              vjust = 0.65,
              size = 0.6,
              textsize = 10,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 0.5, 
                  xmin = 0.85, 
                  xmax = 3.15,
                  annotations = "***"), 
              tip_length = 0.02, 
              manual = T,
              vjust = 0.65,
              size = 0.6,
              textsize = 10,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), 
                     limits = c(0, 0.5), 
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  labs(x = "",
       y = expression("Testicular T (ng mg"^{-1}*")")) +
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
    legend.position = "none")

# Save plot as SVG file
ggsave(filename = "main_figures/plots/Testicular_Testo.svg", plot = p1, 
       width = 6.5, height = 9, bg = "white", dpi = "print", device = "svg")

# Create ggplot object
p2 <- data %>%
  ggplot(aes(fill = morph, y = gon_A4, x = morph)) + 
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
             size = 5, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  geom_signif(data = data,
              aes(y_position = 0.1325, 
                  xmin = 0.85, 
                  xmax = 2.15,
                  annotations = "***"), 
              tip_length = 0.02, 
              manual = T,
              vjust = 0.65,
              size = 0.6,
              textsize = 10,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 0.14, 
                  xmin = 0.85, 
                  xmax = 3.15,
                  annotations = "***"), 
              tip_length = 0.02, 
              manual = T,
              vjust = 0.65,
              size = 0.6,
              textsize = 10,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent","Satellite","Faeder"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent","Satellite","Faeder"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14), 
                     limits = c(0, 0.14), 
                     labels = c("0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  labs(x = "",
       y = expression("Testicular A4 (ng mg"^{-1}*")")) +
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
    legend.position = "none")

# Save plot as SVG file
ggsave(filename = "main_figures/plots/Testicular_A4.svg", plot = p2, 
       width = 6.5, height = 9, bg = "white", dpi = "print", device = "svg")
