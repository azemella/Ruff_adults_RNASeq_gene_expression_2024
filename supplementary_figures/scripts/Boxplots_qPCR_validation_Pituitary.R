# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory 
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("gridExtra")
library("ggplot2")
library("viridis")
library("dplyr")
library("grid")

# Read normalized gene expression count data (RNA-Seq)
PIT_counted <- read.csv(file = "DESeq2_gene_expression_analysis/normalized_counts/normalized_counts_Pit.csv", header = TRUE, row.names = 1)

# Transpose data matrices and convert to data frames
PIT_counted <- as.data.frame(t(PIT_counted))

# Create Morph column for both data frames
PIT_counted <- PIT_counted %>%
  mutate(Morph = case_when(
    grepl("_F_", rownames(PIT_counted)) ~ "Fae",
    grepl("_S_", rownames(PIT_counted)) ~ "Sat",
    grepl("_I_", rownames(PIT_counted)) ~ "Ind",
    TRUE ~ NA_character_
  ))

# Read normalized gene expression count data (qPCR)
qPCR_PIT <- read.csv("data/quantitativePCR/qPCR_pituitary_expression.csv")

# Create ggplot object for PGR expression (RNA-Seq)
p1 <- PIT_counted %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = PGR, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(1000, 2000, 3000, 4000),
                     labels = c("1000", "2000", "3000", "4000"),
                     limits = c(800, 4000)) +
  labs(y = expression(italic("PGR") ~ "expression (RNA-Seq)")) +
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
ggsave(filename = "supplementary_figures/plots/PGR_gene_expression_Pituitary_RNASeq.svg", plot = p1, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")

# Create ggplot object for FSHB expression (RNA-Seq)
p2 <- PIT_counted %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = FSHB, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000),
                     labels = c("0", "5000", "10000", "15000", "20000", "25000"),
                     limits = c(0, 25000)) +
  labs(y = expression(italic("FSHB") ~ "expression (RNA-Seq)")) +
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
ggsave(filename = "supplementary_figures/plots/FSHB_gene_expression_Pituitary_RNASeq.svg", plot = p2, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")

# Create ggplot object for LOC106894170 (GnRHR3) expression (RNA-Seq)
p3 <- PIT_counted %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = LOC106894170, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000),
                     labels = c("0", "1000", "2000", "3000"),
                     limits = c(0, 3500)) +
  labs(y = expression(italic("GnRHR3") ~ "expression (RNA-Seq)")) +
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
ggsave(filename = "supplementary_figures/plots/GnRHR3_gene_expression_Pituitary_RNASeq.svg", plot = p3, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")

# Create ggplot object for PGR expression (qPCR)
p4 <- qPCR_PIT %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = PGR, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6"),
                     limits = c(0, 0.6)) +
  labs(y = expression(italic("PGR") ~ "expression (qPCR)")) +
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
ggsave(filename = "supplementary_figures/plots/PGR_gene_expression_Pituitary_qPCR.svg", plot = p4, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")

# Create ggplot object for FSHB expression (qPCR)
p5 <- qPCR_PIT %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = FSHB, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20),
                     labels = c("0", "5", "10", "15", "20"),
                     limits = c(0, 20)) +
  labs(y = expression(italic("FSHB") ~ "expression (qPCR)")) +
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
ggsave(filename = "supplementary_figures/plots/FSHB_gene_expression_Pituitary_qPCR.svg", plot = p5, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")

# Create ggplot object for GnRHR3 expression (qPCR)
p6 <- qPCR_PIT %>%
  mutate(Morph = factor(Morph, levels = c("Ind", "Sat", "Fae"))) %>%
  ggplot(aes(y = GnRHR3, x = Morph)) + 
  geom_point(size = -1, aes(fill = Morph, color = Morph)) +
  geom_boxplot(aes(fill = Morph),
               position = position_dodge(0.9), 
               width = 0.8, 
               size = 0.35,
               color = "black",
               linewidth = 0.7,
               outlier.shape = NA,
               show.legend = FALSE) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill = Morph), 
               outlier.shape = NA,
               size = 0.35,
               color = "black",
               width = 0.8,
               linewidth = 0.7,
               position = position_dodge(0.9),
               show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., fill = Morph), linewidth = 0.7, width = 0.35,
               size = 0.35, alpha = 1, position = position_dodge(0.9), outlier.shape = NA) +
  geom_point(aes(fill = Morph), 
             size = 6, 
             shape = 21, 
             position = position_jitterdodge(jitter.width = 0.5, 
                                             dodge.width = 0.9),
             stroke = 1) +
  scale_fill_manual(labels = c("Ind", "Sat", "Fae"),
                    values = c("blue", "violetred", "orange")) +
  scale_color_manual(labels = c("Ind", "Sat", "Fae"),
                     values = c("blue", "violetred", "orange")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4),
                     labels = c("0", "0.1", "0.2", "0.3", "0.4"),
                     limits = c(0, 0.4)) +
  labs(y = expression(italic("GnRHR3") ~ "expression (qPCR)")) +
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
ggsave(filename = "supplementary_figures/plots/GnRHR3_gene_expression_Pituitary_qPCR.svg", plot = p6, 
       width = 7, height = 8, bg = "white", dpi = "print", device = "svg")
