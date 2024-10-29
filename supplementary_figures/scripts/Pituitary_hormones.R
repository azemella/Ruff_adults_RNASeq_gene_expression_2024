# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load necessary libraries
library("dplyr")
library("gridExtra")
library("grid")
library("ggplot2")
library("ggpubr")
library("svglite")
library("scales")
library("ggsignif")

# Read normalized count data for pituitary
pit_count_data <- read.csv("DESeq2_gene_expression_analysis/normalized_counts/normalized_counts_Pit.csv", row.names = 1)

# Transpose data frame
pit_count_data_t <- as.data.frame(t(pit_count_data))

# Add phenotype column
pit_count_data_t <- pit_count_data_t %>%
  mutate(Morph = case_when(
    grepl("_F_", names(pit_count_data)) ~ "Faeder",
    grepl("_S_", names(pit_count_data)) ~ "Satellite",
    grepl("_I_", names(pit_count_data)) ~ "Independent",
    TRUE ~ NA_character_
  ))

p1 <- pit_count_data_t %>%
  mutate(Morph = factor(Morph, levels = c("Independent", "Satellite", "Faeder"))) %>%
  ggplot(aes(fill = Morph, y = FSHB, x = Morph)) + 
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
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"), 
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000), 
                     labels = c("0", "5000", "10000", "15000", "20000", "25000"),
                     limits = c(0, 25000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  ggtitle("") +
  xlab("") +
  ylab("FSHB normalized gene expression") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial"),
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
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines")
  ) 

p2 <- pit_count_data_t %>%
  mutate(Morph = factor(Morph, levels = c("Independent", "Satellite", "Faeder"))) %>%
  ggplot(aes(fill = Morph, y = TSHB, x = Morph)) + 
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
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"), 
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 10000, 20000, 30000, 40000), 
                     labels = c("0", "10000", "20000", "30000", "40000"),
                     limits = c(0, 40000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  ggtitle("") +
  xlab("") +
  ylab("TSHB normalized gene expression") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial"),
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
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines")
  ) 

p3 <- pit_count_data_t %>%
  mutate(Morph = factor(Morph, levels = c("Independent", "Satellite", "Faeder"))) %>%
  ggplot(aes(fill = Morph, y = PRL, x = Morph)) + 
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
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"), 
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 400000, 800000, 1200000, 1600000), 
                     labels = c("0", "400000", "800000", "1200000", "1600000"),
                     limits = c(0, 1600000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  ggtitle("") +
  xlab("") +
  ylab("PRL normalized gene expression") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial"),
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
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines")
  ) 

p4 <- pit_count_data_t %>%
  mutate(Morph = factor(Morph, levels = c("Independent", "Satellite", "Faeder"))) %>%
  ggplot(aes(fill = Morph, y = GH1, x = Morph)) + 
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
  geom_signif(data = pit_count_data_t,
              aes(y_position = 4000000, 
                  xmin = 0.75, 
                  xmax = 3.25,
                  annotations = "***"), 
              tip_length = 0.01, 
              manual = T,
              vjust = 0.65,
              size = 0.35,
              textsize = 7,
              inherit.aes = TRUE) +
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"), 
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 1000000, 2000000, 3000000, 4000000), 
                     labels = c("0", "1000000", "2000000", "3000000", "4000000"),
                     limits = c(0, 4000000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  ggtitle("") +
  xlab("") +
  ylab("GH1 normalized gene expression") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial"),
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
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines")
  ) 

p5 <- pit_count_data_t %>%
  mutate(Morph = factor(Morph, levels = c("Independent", "Satellite", "Faeder"))) %>%
  ggplot(aes(fill = Morph, y = POMC, x = Morph)) + 
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
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Independent", "Satellite", "Faeder"), 
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 500000, 1000000, 1500000), 
                     labels = c("0", "500000", "1000000", "1500000"),
                     limits = c(0, 1500000)) +
  guides(color = guide_legend(override.aes = list(shape = 21, size = 4))) +
  ggtitle("") +
  xlab("") +
  ylab("POMC normalized gene expression") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black", family = "Arial"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18, family = "Arial"),
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
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.spacing = unit(0.3, "lines"),
    panel.spacing = unit(0.3, "lines")
  ) 

# Save plot as SVG
svglite("supplementary_figures/plots/Pituitary_hormones.svg", 
        width = 12.5, height = 5)

ggarrange(p1, p2, p3, p4, p5, nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")

dev.off()
