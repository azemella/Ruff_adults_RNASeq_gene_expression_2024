# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory to where the data is located
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("grid")
library("gridExtra")
library("MASS")
library("lattice")
library("latticeExtra")
library("png")
library("RColorBrewer")
library("reshape")
library("reshape2")
library("dplyr")
library("svglite")

# Read ASE results for testes of faeder morph
data_gon_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Gon_ASE.csv")
data_gon_fae <- data_gon_fae[c(1,2,6)]
data_gon_fae$inv_allele_ratio <- data_gon_fae$inv_allele_ratio * 100
data_gon_fae <- data_gon_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for testes of satellite morph
data_gon_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Gon_ASE.csv")
data_gon_sat <- data_gon_sat[c(1,2,6)]
data_gon_sat$inv_allele_ratio <- data_gon_sat$inv_allele_ratio * 100
data_gon_sat <- data_gon_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for hypothalamus tissues of faeder morph
data_hyp_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Hyp_ASE.csv")
data_hyp_fae <- data_hyp_fae[c(1,2,6)]
data_hyp_fae$inv_allele_ratio <- data_hyp_fae$inv_allele_ratio * 100
data_hyp_fae <- data_hyp_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for hypothalamus of satellite morph
data_hyp_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Hyp_ASE.csv")
data_hyp_sat <- data_hyp_sat[c(1,2,6)]
data_hyp_sat$inv_allele_ratio <- data_hyp_sat$inv_allele_ratio * 100
data_hyp_sat <- data_hyp_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for lateral septum of faeder morph
data_ls_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_LS_ASE.csv")
data_ls_fae <- data_ls_fae[c(1,2,6)]
data_ls_fae$inv_allele_ratio <- data_ls_fae$inv_allele_ratio * 100
data_ls_fae <- data_ls_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for lateral septum of satellite morph
data_ls_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_LS_ASE.csv")
data_ls_sat <- data_ls_sat[c(1,2,6)]
data_ls_sat$inv_allele_ratio <- data_ls_sat$inv_allele_ratio * 100
data_ls_sat <- data_ls_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for caudal nidopallium of faeder morph
data_nc_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_NC_ASE.csv")
data_nc_fae <- data_nc_fae[c(1,2,6)]
data_nc_fae$inv_allele_ratio <- data_nc_fae$inv_allele_ratio * 100
data_nc_fae <- data_nc_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for caudal nidopallium of satellite morph
data_nc_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_NC_ASE.csv")
data_nc_sat <- data_nc_sat[c(1,2,6)]
data_nc_sat$inv_allele_ratio <- data_nc_sat$inv_allele_ratio * 100
data_nc_sat <- data_nc_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for pituitary gland of faeder morph
data_pit_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Pit_ASE.csv")
data_pit_fae <- data_pit_fae[c(1,2,6)]
data_pit_fae$inv_allele_ratio <- data_pit_fae$inv_allele_ratio * 100
data_pit_fae <- data_pit_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for pituitary gland of satellite morph
data_pit_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Pit_ASE.csv")
data_pit_sat <- data_pit_sat[c(1,2,6)]
data_pit_sat$inv_allele_ratio <- data_pit_sat$inv_allele_ratio * 100
data_pit_sat <- data_pit_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for preoptic area of faeder morph
data_poa_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_POA_ASE.csv")
data_poa_fae <- data_poa_fae[c(1,2,6)]
data_poa_fae$inv_allele_ratio <- data_poa_fae$inv_allele_ratio * 100
data_poa_fae <- data_poa_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for preoptic area of satellite morph
data_poa_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_POA_ASE.csv")
data_poa_sat <- data_poa_sat[c(1,2,6)]
data_poa_sat$inv_allele_ratio <- data_poa_sat$inv_allele_ratio * 100
data_poa_sat <- data_poa_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for raphe of faeder morph
data_rap_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Rap_ASE.csv")
data_rap_fae <- data_rap_fae[c(1,2,6)]
data_rap_fae$inv_allele_ratio <- data_rap_fae$inv_allele_ratio * 100
data_rap_fae <- data_rap_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for raphe of satellite morph
data_rap_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Rap_ASE.csv")
data_rap_sat <- data_rap_sat[c(1,2,6)]
data_rap_sat$inv_allele_ratio <- data_rap_sat$inv_allele_ratio * 100
data_rap_sat <- data_rap_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for nucleus taeniae and arcopallium of faeder morph
data_tna_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_A+TnA_ASE.csv")
data_tna_fae <- data_tna_fae[c(1,2,6)]
data_tna_fae$inv_allele_ratio <- data_tna_fae$inv_allele_ratio * 100
data_tna_fae <- data_tna_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for nucleus taeniae and arcopallium of satellite morph
data_tna_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_A+TnA_ASE.csv")
data_tna_sat <- data_tna_sat[c(1,2,6)]
data_tna_sat$inv_allele_ratio <- data_tna_sat$inv_allele_ratio * 100
data_tna_sat <- data_tna_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for ventral tegmental area and substantia nigra of faeder morph
data_vta_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_VTA+SN_ASE.csv")
data_vta_fae <- data_vta_fae[c(1,2,6)]
data_vta_fae$inv_allele_ratio <- data_vta_fae$inv_allele_ratio * 100
data_vta_fae <- data_vta_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for ventral tegmental area and substantia nigra of satellite morph
data_vta_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_VTA+SN_ASE.csv")
data_vta_sat <- data_vta_sat[c(1,2,6)]
data_vta_sat$inv_allele_ratio <- data_vta_sat$inv_allele_ratio * 100
data_vta_sat <- data_vta_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for adrenal glands of faeder morph
data_adr_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Adr_ASE.csv")
data_adr_fae <- data_adr_fae[c(1,2,6)]
data_adr_fae$inv_allele_ratio <- data_adr_fae$inv_allele_ratio * 100
data_adr_fae <- data_adr_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for adrenal glands of satellite morph
data_adr_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Adr_ASE.csv")
data_adr_sat <- data_adr_sat[c(1,2,6)]
data_adr_sat$inv_allele_ratio <- data_adr_sat$inv_allele_ratio * 100
data_adr_sat <- data_adr_sat %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for liver of faeder morph
data_liv_fae <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/faeder/results_tissue_faeder_Liv_ASE.csv")
data_liv_fae <- data_liv_fae[c(1,2,6)]
data_liv_fae$inv_allele_ratio <- data_liv_fae$inv_allele_ratio * 100
data_liv_fae <- data_liv_fae %>% mutate_if(is.numeric, ~round(., 1))

# Read ASE results for liver of satellite morph
data_liv_sat <- read.csv(file = "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Liv_ASE.csv")
data_liv_sat <- data_liv_sat[c(1,2,6)]
data_liv_sat$inv_allele_ratio <- data_liv_sat$inv_allele_ratio * 100
data_liv_sat <- data_liv_sat %>% mutate_if(is.numeric, ~round(., 1))

# Combine data for all tissues for each morph separately
data_sat <- bind_rows(data_gon_sat, data_hyp_sat, data_ls_sat, data_nc_sat, data_pit_sat, data_poa_sat, data_rap_sat, data_tna_sat, data_vta_sat, data_liv_sat, data_adr_sat)
data_fae <- bind_rows(data_gon_fae, data_hyp_fae, data_ls_fae, data_nc_fae, data_pit_fae, data_poa_fae, data_rap_fae, data_tna_fae, data_vta_fae, data_liv_fae, data_adr_fae)

# Convert data to matrices for plotting heatmaps
mat_cast_sat <- dcast(data_sat, gene ~ Tissue, value.var = "inv_allele_ratio")
data_matrix_sat <- as.matrix(mat_cast_sat[, -c(1)])
rownames(data_matrix_sat) <- mat_cast_sat$gene

mat_cast_fae <- dcast(data_fae, gene ~ Tissue, value.var = "inv_allele_ratio")
data_matrix_fae <- as.matrix(mat_cast_fae[, -c(1)])
rownames(data_matrix_fae) <- mat_cast_fae$gene

# Rename and reorder columns for faeder morph
colnames(data_matrix_fae) <- c("ADR", "NC", "TES", "HYP", "LS", "LIV", "A+TnA", "PIT", "POM", "RAP", "VTA+SN")
data_matrix_fae <- data_matrix_fae[, c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV")]

# Rename and reorder columns for satellite morph
colnames(data_matrix_sat) <- c("ADR", "NC", "TES", "HYP", "LS", "LIV", "A+TnA", "PIT", "POM", "RAP", "VTA+SN")
data_matrix_sat <- data_matrix_sat[, c("A+TnA", "LS", "POM", "HYP", "NC", "RAP", "VTA+SN", "PIT", "TES", "ADR", "LIV")]

# Remove specific genes (Remove CENPN and non-coding RNAs) from both data matrices
row_to_remove_fae <- which(row.names(data_matrix_fae) %in% c("CENPN", "LOC106888090", "LOC106888214", "LOC106888201", "LOC106888296"))
row_to_remove_sat <- which(row.names(data_matrix_sat) %in% c("CENPN", "LOC106888090", "LOC106888214", "LOC106888201", "LOC106888296"))
data_matrix_fae <- data_matrix_fae[-row_to_remove_fae, ]
data_matrix_sat <- data_matrix_sat[-row_to_remove_sat, ]

# Define genes that are not expressed in all tissues to be added
new_genes_fae <- c("MAP1LC3B", "TUBB3", "DEF8")
new_genes_sat <- c("ZDHHC7", "LOC106888097", "FOXL1", "FOXC2", "FOXF1", "IRF8", "COX4I1", "FBXO31", "MAP1LC3B", "ZNF469", "ZFPM1", "CYBA", 
                   "RNF166", "CTU2", "APRT", "GALNS", "TRAPPC2L", "CBFA2T3", "DEF8", "DBNDD1")

# Add new genes to the data matrices
data_matrix_fae <- rbind(data_matrix_fae, matrix(NA, nrow = length(new_genes_fae), ncol = ncol(data_matrix_fae)))
rownames(data_matrix_fae)[(nrow(data_matrix_fae) - length(new_genes_fae) + 1):nrow(data_matrix_fae)] <- new_genes_fae

data_matrix_sat <- rbind(data_matrix_sat, matrix(NA, nrow = length(new_genes_sat), ncol = ncol(data_matrix_sat)))
rownames(data_matrix_sat)[(nrow(data_matrix_sat) - length(new_genes_sat) + 1):nrow(data_matrix_sat)] <- new_genes_sat

# Rename a some genes
old_row_names <- c("LOC106888173", "DUSP15", "LOC106888312", "LOC106888185", "LOC106888094", "LOC106888384")
new_row_names <- c("CDH3", "DUSP22A-like", "HNF4B", "SULT2B1", "CDH1-like", "CDH1")

rows_to_replace <- match(old_row_names, rownames(data_matrix_sat))
rownames(data_matrix_sat)[rows_to_replace] <- new_row_names

rows_to_replace <- match(old_row_names, rownames(data_matrix_fae))
rownames(data_matrix_fae)[rows_to_replace] <- new_row_names

# Read inversion genes list
inv_genes <- read.csv(file = "metadata/NCBI_inversion_genes_list.csv")

# Convert to data frame and set gene names as row names
inv_genes <- as.data.frame(inv_genes$Gene_symbol)
colnames(inv_genes)[1] <- "Gene_ID"
rownames(inv_genes) <- inv_genes$Gene_ID

rows_to_remove <- c("CENPN", "LOC106888090", "LOC106888214", "LOC106888201", "LOC106888296")
inv_genes <- as.data.frame(inv_genes[!(rownames(inv_genes) %in% rows_to_remove), ])
colnames(inv_genes)[1] <- "Gene_ID"
rownames(inv_genes) <- inv_genes$Gene_ID

# Sort data matrices
sorted_data_matrix_fae <- data_matrix_fae[rownames(inv_genes), colnames(data_matrix_fae)]
sorted_data_matrix_sat <- data_matrix_sat[rownames(inv_genes), colnames(data_matrix_sat)]

# Flip data matrices' inversion genes order
sorted_data_matrix_fae <- sorted_data_matrix_fae[rev(rownames(sorted_data_matrix_fae)), ]
sorted_data_matrix_sat <- sorted_data_matrix_sat[rev(rownames(sorted_data_matrix_sat)), ]

# Transpose data matrices
sorted_data_matrix_fae <- t(sorted_data_matrix_fae)
sorted_data_matrix_sat <- t(sorted_data_matrix_sat)

# Define colors for heatmap
colours <- grDevices::colorRampPalette(
  colors = c('#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B', '#DD3D2D', '#A50026'))(100)

# Define numbers for axis ticks
numbers_h <- seq(1.5, 120, by = 1)
numbers_v <- seq(1.5, 120, by = 1)

# Set padding for axis ticks
lattice.options(axis.padding = list(factor = 0.5))

# Set up SVG output file for faeder morph heatmap
svglite(filename = "supplementary_figures/plots/Heatmap_ASE_inv_morphs_11tissues_Fae.svg", 
        width = 8,
        height = 20, 
        pointsize = 16,
        bg = "white")

# Create heatmap for the faeder morph
levelplot(sorted_data_matrix_fae,
          xlab = NULL,
          ylab = NULL,
          col.regions = colours,
          margin = TRUE,
          at = seq.int(from = 0, to = 100, by = 5),
          label.style = "flat",
          colorkey = list(at = seq(0, 100, by = 5),
                          col = colours,
                          space = "left", height = 0.5, width = 1.5,
                          labels = list(cex = 1.75,
                                        font = 1,
                                        at = c(0, 20, 40, 60, 80, 100), 
                                        labels = c("0", "20", "40", "60", "80", "100"))),
          aspect = "fill",
          par.settings = list(axis.text = list(fontfamily = "arial"),
                              par.xlab.text = list(fontfamily = "arial"),
                              par.ylab.text = list(fontfamily = "arial"),
                              layout.widths = list(left.padding = 3)),
          scales = list(x = list(cex = 2.25, alternating = 2, font = 1, tck = c(0,1), rot = 90), 
                        y = list(cex = 0.6, alternating = 0, tck = c(0,1))),
          panel = function(...){
            panel.levelplot(...)
            panel.abline(v = c(numbers_v), lwd = 1.5)
            panel.abline(h = c(numbers_h), lwd = 0.25)})

# Close SVG device
dev.off()

# Set up SVG output file for satellite morph heatmap
svglite(filename = "supplementary_figures/plots/Heatmap_ASE_inv_morphs_11tissues_Sat.svg", 
        width = 8,
        height = 20, 
        pointsize = 16,
        bg = "white")

# Create heatmap for the satellite morph
levelplot(sorted_data_matrix_sat,
          xlab = NULL,
          ylab = NULL,
          col.regions = colours,
          margin = TRUE,
          at = seq.int(from = 0, to = 100, by = 5),
          label.style = "flat",
          colorkey = list(at = seq(0, 100, by = 5),
                          col = colours,
                          space = "right", height = 0.5, width = 1.5,
                          labels = list(cex = 1.75,
                                        font = 1,
                                        at = c(0, 20, 40, 60, 80, 100), 
                                        labels = c("0", "20", "40", "60", "80", "100"))),
          aspect = "fill",
          par.settings = list(axis.text = list(fontfamily = "arial"),
                              par.xlab.text = list(fontfamily = "arial"),
                              par.ylab.text = list(fontfamily = "arial"),
                              layout.widths = list(left.padding = 3)),
          scales = list(x = list(cex = 2.25, alternating = 2, font = 1, tck = c(0,1), rot = 90), 
                        y = list(cex = 0.6, alternating = 0, tck = c(1,0))),
          panel = function(...){
            panel.levelplot(...)
            panel.abline(v = c(numbers_v), lwd = 1.5)
            panel.abline(h = c(numbers_h), lwd = 0.25)})

# Close SVG device
dev.off()
