# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("reshape2")

# Read kinetics data
melted_data <- read.csv(file = "data/kinetics/Data_MM.csv")

# Convert 'Morph' column to a factor with specified levels
melted_data$Morph <- factor(melted_data$Morph, levels = c("IND", "SAT", "FAE", "VAR1", "VAR2", "VAR3", "VAR4"))

# Define colors, fills, and shapes for plotting
cols <- c("blue", "violetred", "orange", 'darkgrey', 'darkgrey', 'darkgrey', 'darkgrey')
fills <- c("blue", "violetred", "orange", 'darkgrey', 'darkgrey', 'darkgrey', 'darkgrey')
shapes <- c(19, 19, 19, 19, 24, 25, 22)

# Create ggplot object
p <- ggplot(melted_data, aes(x = Substrate, y = Concentration)) +
  labs(x = expression("Testosterone concentration (nM)"),
       y = expression("Testosterone conversion rate (pmol min"^{-1}*")  ")) +
  theme_bw() +
  geom_point(data = subset(melted_data, Morph %in% c("VAR1", "VAR2", "VAR3", "VAR4")),
             aes(color = Morph, shape = Morph, fill = Morph), 
             size = 3.5) +
  geom_point(data = subset(melted_data, Morph %in% c("IND", "SAT", "FAE")),
             aes(color = Morph, shape = Morph, fill = Morph), 
             size = 3.5) +
  geom_errorbar(aes(ymin = Concentration-STD, ymax = Concentration+STD, color = Morph, shape = Morph, fill = Morph), 
                width=.2,
                alpha = 0.5,
                position = position_dodge(0.05)) +
  scale_color_manual(values = c("IND" = cols[1], 'SAT' = cols[2], 'FAE' = cols[3],
                                "VAR1" = cols[4], "VAR2" = cols[5], "VAR3" = cols[6], "VAR4" = cols[7]),
                     breaks = c("IND", "SAT", "FAE", "VAR1", "VAR2", "VAR3", "VAR4"),
                     labels = c("IND", "SAT", "FAE", "VAR-1", "VAR-2", "VAR-3", "VAR-4")) +
  scale_shape_manual(values = c("IND" = shapes[1], 'SAT' = shapes[2], 'FAE' = shapes[3],
                                "VAR1" = shapes[4], "VAR2" = shapes[5], "VAR3" = shapes[6], "VAR4" = shapes[7]),
                     breaks = c("IND", "SAT", "FAE", "VAR1", "VAR2", "VAR3", "VAR4"),
                     labels = c("IND", "SAT", "FAE", "VAR-1", "VAR-2", "VAR-3", "VAR-4")) +
  scale_fill_manual(values = c("IND" = fills[1], 'SAT' = fills[2], 'FAE' = fills[3],
                               "VAR1" = fills[4], "VAR2" = fills[5], "VAR3" = fills[6], "VAR4" = fills[7]),
                    breaks = c("IND", "SAT", "FAE", "VAR1", "VAR2", "VAR3", "VAR4"),
                    labels = c("IND", "SAT", "FAE", "VAR-1", "VAR-2", "VAR-3", "VAR-4")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), 
                     limits = c(0, 3), 
                     labels = c("0", "1", "2", "3")) +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[4], data = filter(melted_data, Morph == "VAR1"), linetype="dashed") +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[5], data = filter(melted_data, Morph == "VAR2"), linetype="dashed") +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[6], data = filter(melted_data, Morph == "VAR3"), linetype="dashed") +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[7], data = filter(melted_data, Morph == "VAR4"), linetype="dashed") +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[1], data = filter(melted_data, Morph == "IND")) +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[2], data = filter(melted_data, Morph == "SAT")) +
  geom_smooth(method = "nls", formula = y ~ Vmax * x / (Km + x), start = list(Vmax = 50, Km = 2),
              se = F, size = 1.25, color = cols[3], data = filter(melted_data, Morph == "FAE")) +
  theme(axis.title.x.bottom = element_text(size = 26, family = "Arial",margin = margin(15, 0, 0, 0)),
        axis.title.y.left = element_text(size = 26, family = "Arial",margin = margin(0, 15, 0, 0)),
        axis.text.x.bottom = element_text(size = 18, colour = "black", family = "Arial"),
        axis.text.y.left = element_text(size = 18, colour = "black", family = "Arial"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.75),
        axis.ticks.x = element_line(color = "black", linewidth = 0.75),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = alpha("white", 1)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 24, colour = "black", family = "Arial")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5, order = c("IND", "SAT", "FAE", "VAR1", "VAR2", "VAR3", "VAR4"))))

# Save plot as SVG file
ggsave(filename = "main_figures/plots/MM_plot.svg", plot = p, 
       width = 10, height = 7, bg = "white", dpi = "print", device = "svg")
