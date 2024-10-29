# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

### Open both the two R scripts "Hormone_plots.R" and "Behavior_plots.R" on Rstudio.  
### Run "Hormone_plots.R" first, followed by "Behavior_plots.R". The output is a single SVG file containing
### a total of four boxplots: two showing blood hormone measurements and two showing behavioral quantifications

# Set the working directory to where the data files are located
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("dplyr")
library("gridExtra")
library("grid")
library("ggplot2")
library("ggpubr")
library("svglite")

# Read data files for male and female ruffs
data <- read.csv("data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv", header = T)

# Remove unnecessary columns
data <- data[, -c(5:12)]

# Rename columns
colnames(data)[2] <- "Morph"
colnames(data)[3] <- "Testo"
colnames(data)[4] <- "A4"

# Convert hormone values to ng/ml
data$Testo <- data$Testo / 1000
data$A4 <- data$A4 / 1000
data$Morph <- factor(data$Morph, levels = c("Ind", "Sat", "Fae"))

# Generate plot for circulating testosterone levels
plotTT <- data %>%
  ggplot(aes(fill = Morph, y = Testo, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.5,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               size = 0.35,
               width = 0.8,
               linewidth = 0.5,
               position = position_dodge(0.9),
               show.legend = FALSE,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = data,
              aes(y_position = 8.25, 
                  xmin = 0.85, 
                  xmax = 2.25,
                  annotations = "**"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 8.75, 
                  xmin = 0.85, 
                  xmax = 3.2,
                  annotations = "**"), 
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
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), 
                     limits = c(0, 8.75), 
                     labels = c("0", "2", "4", "6", "8")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  labs(x = "",
       y = expression("Circulating T (ng ml"^{-1}*")")) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11.5, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial", margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white", 1)),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box.just = "left",
    legend.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 20, family = "Arial"))

# Generate plot for circulating A4 levels
plotA4 <- data %>%
  ggplot(aes(fill = Morph, y = A4, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph), 
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               linewidth = 0.5,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               size = 0.35,
               width = 0.8,
               linewidth = 0.5,
               position = position_dodge(0.9),
               show.legend = FALSE,
               outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.5, width = 0.5, size = 0.35, alpha = 1, 
               position = position_dodge(0.9), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 2, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9)) +
  geom_signif(data = data,
              aes(y_position = 23.5, 
                  xmin = 0.85, 
                  xmax = 2.25,
                  annotations = "***"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  geom_signif(data = data,
              aes(y_position = 25, 
                  xmin = 0.85, 
                  xmax = 3.2,
                  annotations = "***"), 
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
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25), 
                     limits = c(0, 25), 
                     labels = c("0", "5", "10", "15", "20", "25")) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 6))) +
  labs(x = "",
       y = expression("Circulating A4 (ng ml"^{-1}*")")) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 11.5, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial", margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = alpha("white")),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, family = "Arial"),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box.just = "left",
    legend.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "white"),
    strip.text.x = element_text(size = 20, family = "Arial"))
