# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("ggplot2")
library("dplyr")

# Read data frame with number of DE inversion genes
data <- read.csv("others/DESeq2_total_number_of_DEGs.csv")
data <- data[, -3]

# Define custom order for tissues and comparisons
custom_order <- c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV")
custom_order2 <- c("Fae vs. Ind", "Sat vs. Ind", "Fae vs. Sat")

# Apply custom order to the factors
data$Tissue <- factor(data$Tissue, levels = custom_order)
data$Comparison <- factor(data$Comparison, levels = custom_order2)

# Create an empty data frame for plotting points
df <- data.frame(
  areas = as.numeric(c(NA, NA, NA)),
  Comparison = c("Fae vs. Ind", "Sat vs. Ind", "Fae vs. Sat"),
  Value = as.numeric(c(NA, NA, NA))
)

# Create ggplot object
p <- ggplot(data, aes(x = Tissue, y = Inversion_DEGs, fill = Comparison)) +
  geom_abline(slope = 0, intercept = 5,  col = "darkgrey", lty = 2) +
  geom_abline(slope = 0, intercept = 10,  col = "darkgrey", lty = 2) +
  geom_abline(slope = 0, intercept = 15,  col = "darkgrey", lty = 2) +
  geom_abline(slope = 0, intercept = 20,  col = "darkgrey", lty = 2) +
  geom_bar(stat = "identity", 
           color = "black", 
           position = position_dodge(width = 0.75), 
           width = 0.75,
           show.legend = FALSE) +
  labs(title = "", x = "", y = "Number of DEGs in\n inverted regions") +
  scale_fill_manual(labels = c("Fae vs. Ind", "Sat vs. Ind", "Fae vs. Sat"), 
                    values = c("black", "grey", "white")) +
  scale_color_manual(labels = c("Fae vs. Ind", "Sat vs. Ind", "Fae vs. Sat"), 
                     values = c("black", "grey", "white")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 21)) +
  geom_point(data = df, aes(x = areas, y = Value, fill = Comparison), color = "black", size = 9, shape = 21, stroke = 0.75) +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", size = 28, family = "Arial", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(color = "black", size = 18, family = "Arial"),
        axis.ticks.x.bottom = element_line(linewidth = 0.75, color = "black"),
        axis.ticks.y.left = element_line(linewidth = 0.75, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.x.bottom = element_line(color = 'black', linewidth = 0.75),
        axis.line.y.left   = element_line(color = 'black', linewidth = 0.75),
        axis.line.y.right  = element_blank(),
        axis.text.y.right  = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 30, family = "Arial", margin = margin(r = 10)),
        legend.position = "top",
        axis.title.y = element_text(color = "black", size = 30, family = "Arial", hjust = 0.65, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.75, "cm")) 

# Save plot as SVG file
ggsave(filename = "supplementary_figures/plots/Barplot_number_inversion_genes.svg", plot = p, 
       width = 10, height = 8, bg = "white", dpi = "print", device = "svg")
