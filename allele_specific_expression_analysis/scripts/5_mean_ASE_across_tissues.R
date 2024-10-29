# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory and define directory path
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Run this script for "faeder" and "satellite" files separately
# Replace the inversion morph name in lines #12, #16, #20, #24, #28, #32, #36, #40, #44, #48, #52, #56 and #67

# Read allele-specific expression (ASE) data for each tissue and rename columns
ase_tna <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_A+TnA_ASE.csv"))
ase_tna <- ase_tna[,c(1,6)]
colnames(ase_tna)[2] <- "inv_allele_ratio_tna"

ase_adr <- read.csv(file.path(dir, "allele_specific_expression/ASE_results/tissues/satellite/results_tissue_satellite_Adr_ASE.csv"))
ase_adr <- ase_adr[,c(1,6)]
colnames(ase_adr)[2] <- "inv_allele_ratio_adr"

ase_blo <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Blood_ASE.csv"))
ase_blo <- ase_blo[,c(1,6)]
colnames(ase_blo)[2] <- "inv_allele_ratio_blo"

ase_gon <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Gon_ASE.csv"))
ase_gon <- ase_gon[,c(1,6)]
colnames(ase_gon)[2] <- "inv_allele_ratio_gon"

ase_hyp <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Hyp_ASE.csv"))
ase_hyp <- ase_hyp[,c(1,6)]
colnames(ase_hyp)[2] <- "inv_allele_ratio_hyp"

ase_liv <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Liv_ASE.csv"))
ase_liv <- ase_liv[,c(1,6)]
colnames(ase_liv)[2] <- "inv_allele_ratio_liv"

ase_ls <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_LS_ASE.csv"))
ase_ls <- ase_ls[,c(1,6)]
colnames(ase_ls)[2] <- "inv_allele_ratio_ls"

ase_nc <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_NC_ASE.csv"))
ase_nc <- ase_nc[,c(1,6)]
colnames(ase_nc)[2] <- "inv_allele_ratio_nc"

ase_pit <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Pit_ASE.csv"))
ase_pit <- ase_pit[,c(1,6)]
colnames(ase_pit)[2] <- "inv_allele_ratio_pit"

ase_poa <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_POA_ASE.csv"))
ase_poa <- ase_poa[,c(1,6)]
colnames(ase_poa)[2] <- "inv_allele_ratio_poa"

ase_rap <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_Rap_ASE.csv"))
ase_rap <- ase_rap[,c(1,6)]
colnames(ase_rap)[2] <- "inv_allele_ratio_rap"

ase_vta <- read.csv(file.path(dir, "allele_specific_expression_analysis/ASE_results/tissues/satellite/results_tissue_satellite_VTA+SN_ASE.csv"))
ase_vta <- ase_vta[,c(1,6)]
colnames(ase_vta)[2] <- "inv_allele_ratio_vta"

# Merge data frames and calculate mean allele-specific expression (ASE) for each gene across tissues
df_list <- list(ase_tna, ase_adr, ase_blo, ase_gon, ase_hyp, ase_liv, ase_ls, ase_nc, ase_pit, ase_poa, ase_rap, ase_vta)
mean_ase_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "gene"), df_list)

# Calculate mean ASE and write results to CSV file
mean_ase_df$mean_ase <- rowMeans(mean_ase_df[,2:12], na.rm=TRUE)
write.csv(as.data.frame(mean_ase_df), 
          file = "Mean_ASE_11_tissues_Satellite.csv")
