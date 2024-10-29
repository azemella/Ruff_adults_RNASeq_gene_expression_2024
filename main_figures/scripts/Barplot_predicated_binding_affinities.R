# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load necessary libraries
library("ggplot2")
library("scales")

# Create data frame
data <- data.frame(
  morph = c("Independent", "Independent", "Satellite", "Satellite", "Faeder", "Faeder"),
  ligand = c("Testosterone", "NAD+", "Testosterone", "NAD+", "Testosterone", "NAD+"),
  values = c(8.6, 7.54, 8.99, 7.88, 8.91, 7.94)
)

# Set factor levels for 'morph' and 'ligand' variables
data$morph <- factor(data$morph, levels = c("Independent", "Satellite", "Faeder"))
data$ligand <- factor(data$ligand, levels = c("NAD+", "Testosterone"))

# Create ggplot object
p <- ggplot(data, aes(x = values, 
                       y = ligand, 
                       fill = morph)) +
  geom_vline(xintercept = 6,  col = "darkgrey", linewidth = 0.25, lty = 2) +
  geom_vline(xintercept = 7,  col = "darkgrey", linewidth = 0.25, lty = 2) +
  geom_vline(xintercept = 8,  col = "darkgrey", linewidth = 0.25, lty = 2) +
  geom_vline(xintercept = 9,  col = "darkgrey", linewidth = 0.25, lty = 2) +
  geom_bar(position = "dodge", 
           stat = "identity", 
           aes(fill = morph, color = morph),
           width = 0.7,
           key_glyph = draw_key_point, 
           colour = "black") +
  geom_text(aes(label = values * (-1)), 
            position = position_dodge(width = 0.7), 
            vjust = 0.5,
            hjust = -0.1,
            size = 8) +  
  scale_x_continuous(position = "top", 
                     expand = c(0, 0),
                     breaks = c(5, 6, 7, 8, 9),
                     labels = c("-5", "-6", "-7", "-8", "-9")) +
  scale_y_discrete(labels = c(expression(NAD^"+"), "Testosterone"),
                   expand = c(0.5, 0)) +
  scale_fill_manual(labels = c("Independent", "Satellite", "Faeder"), 
                    values = c("blue", "violetred", "orange")) +
  coord_cartesian(xlim = c(5, 9.5)) +
  labs(title = "", 
       x = expression("Predicted binding affinity \u0394 G Binding (kcal mol-1)"), 
       y = "") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 8, colour = "white"))) +
  theme_bw() +
  theme(
    axis.text.x.top = element_text(color = "black", size = 20, family = "Arial"), 
    axis.text.y.left = element_text(color = "black", size = 30, family = "Arial", angle = 90, hjust = 0.5),
    axis.ticks.x.top = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.y.left = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.title.y = element_blank(),
    axis.title.x.top = element_text(color = "black", size = 30, family = "Arial", hjust = 0.65, margin = margin(t = 0, r = 15, b = 15, l = 15)),
    axis.line.x.top = element_line(color = 'black', linewidth = 0.75),
    axis.line.y.left   = element_line(color = 'black', linewidth = 0.75),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 26, family = "Arial"),
    legend.position = c(.8725, 0.125)
  )

# Save plot as SVG
ggsave(filename = "main_figures/plots/Barplot_binding_affinities.svg", plot = p, 
       width = 13, height = 7, bg = "white", dpi = "print", device = "svg")
