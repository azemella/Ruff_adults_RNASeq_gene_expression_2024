# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("ggplot2")
library("scales")

# Read data from CSV files
res1 <- read.csv(file = "GO_terms_enrichment_analysis/results/clusterProfiler_GO_ORA_faeder_vs_independent_Gon.csv", 
                 header = TRUE)
res2 <- read.csv(file = "GO_terms_enrichment_analysis/results/clusterProfiler_GO_ORA_satellite_vs_independent_Gon.csv", 
                 header = TRUE)
res3 <- read.csv(file = "GO_terms_enrichment_analysis/results/clusterProfiler_GO_ORA_faeder_vs_satellite_Gon.csv", 
                 header = TRUE)

# Select top non-redundant GO terms in at least two lists
res1_keep <- c("GO:1901615", "GO:0034698", "GO:0006590", "GO:0008202", "GO:0008203", "GO:0042445", "GO:0032354", "GO:0010817", "GO:0044283", 
               "GO:0120178", "GO:0033327", "GO:0008584", "GO:0008211")
subset_res1 <- res1[res1$ID %in% res1_keep & res1$ONTOLOGY != "MF", ]

res2_keep <- c("GO:0008211", "GO:0120178", "GO:0034698", "GO:0008584", "GO:0006694", "GO:0042445", "GO:0008202", "GO:0006547", "GO:0033327", 
               "GO:0032354", "GO:1901615", "GO:0010817", "GO:0032350", "GO:0044283")
subset_res2 <- res2[res2$ID %in% res2_keep & res2$ONTOLOGY != "MF", ]

res3_keep <- c("GO:0034728", "GO:0006547", "GO:0051593", "GO:0006338", "GO:0071824")
subset_res3 <- res3[res3$ID %in% res3_keep & res3$ONTOLOGY != "MF", ]

# Add comparison labels to subsets
subset_res1$Comparison <- "Fae vs. Ind"
subset_res2$Comparison <- "Sat vs. Ind"
subset_res3$Comparison <- "Fae vs. Sat"

# Combine subsets into one dataframe
data_df <- rbind(subset_res1, subset_res2, subset_res3)

# Capitalize the first letter of Description column
data_df$Description <- paste0(toupper(substr(data_df$Description, 1, 1)), substring(data_df$Description, 2))

# Create ggplot object
p <- ggplot(data_df, aes(x = factor(data_df$Comparison, levels = c("Fae vs. Sat", "Sat vs. Ind", "Fae vs. Ind")), 
                         y = factor(data_df$Description, levels = c("Response to gonadotropin", "Glucocorticoid metabolic process", "Thyroid hormone generation", 
                                                                    "Steroid hormone biosynthetic process", "Steroid metabolic process", "Hormone metabolic process", 
                                                                    "Cholesterol metabolic process", "Response to follicle-stimulating hormone", "Male gonad development", 
                                                                    "Steroid biosynthetic process", "Regulation of hormone levels", "Small molecule biosynthetic process", 
                                                                    "Nucleosome organization", "Histidine metabolic process", "Leydig cell differentiation", 
                                                                    "Chromatin remodeling", "Organic hydroxy compound metabolic process", "Regulation of hormone metabolic process", 
                                                                    "Protein-DNA complex organization", "Response to folic acid")),
                         fill = FDR, size = Count)) +
  geom_point(aes(size = Count, fill = FDR), colour = "black", pch = 21) +
  scale_fill_gradientn(colors = c("#b3eebe", "#46bac2", "#371ea3"),
                       guide = guide_colorbar(reverse = FALSE, order = 1, ticks.colour = "black", frame.colour = "black",
                                              frame.linewidth = 1/.pt, ticks.linewidth = 1/.pt),
                       limits = c(0.000001, 0.1),
                       breaks = c(0.01, 0.05, 0.095),
                       labels = c("0.01","0.05", "0.1")) +
  scale_size_continuous(range = c(5, 13), 
                        breaks = c(1, 4, 8, 12),
                        name = "Gene count") +
  scale_x_discrete(labels = c("Fae vs. Ind" = "Faeder vs.\n Independent", 
                              "Sat vs. Ind" = "Satellite vs.\n Independent", 
                              "Fae vs. Sat" = "Faeder vs.\n Satellite")) +
  theme_bw() + 
  ylab("") + 
  xlab("") +
  coord_flip() +
  theme(axis.text.x = element_text(size = 16, colour = "black", family = "Arial", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black", family = "Arial", hjust = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.y = unit(.3, "cm"),
        axis.ticks.length.x = unit(.15, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 'dotted', colour = '#808080'),
        panel.background = element_rect(fill = alpha("white", 1)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.title = element_text(size = 16, colour = "black", family = "Arial"),
        legend.text = element_text(size = 14, colour = "black", family = "Arial"),
        legend.title.align = 0.5,
        legend.text.align = 0,
        legend.spacing.x = unit(0.2, "cm"))
  
# Save the plot as SVG file
ggsave(filename = "supplementary_figures/plots/GO_enrichment_gonads.svg", plot = p, 
      width = 12, height = 8, bg = "white", dpi = "print", device = "svg")
