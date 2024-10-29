# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory and define directory path
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")

# Run this script for "faeder" and "satellite" files separately
# Replace the tissue type in lines #17, #60 and #71
# Replace the inversion morph name in lines #16 and #71

# Get a list of files in the specified directory
my_files <- list.files(path = "allele_specific_expression_analysis/Normalization/faeder/", 
                       pattern = "*_HYP_*", full.names = TRUE)

# Read and merge data from each file in the list
for (file in my_files){
  if (!exists("tissue_df")){
    tissue_df <- read.csv(file)
  }
  if (exists("tissue_df")){
    temp_dataset <-read.csv(file)
    tissue_df <- rbind(tissue_df, temp_dataset)
    rm(temp_dataset)
  }
}

# Remove duplicate rows
tissue_df <- tissue_df %>% distinct()

# Summarize counts at the gene level
tissue_df1 <- aggregate(tissue_df$total_count, by = list(gene = tissue_df$gene), FUN = sum)
colnames(tissue_df1)[2] <- "total_count"

tissue_df2 <- aggregate(tissue_df$norm_alt_count, by = list(gene = tissue_df$gene), FUN = sum)
colnames(tissue_df2)[2] <- "alt_allele_count"

tissue_df3 <- aggregate(tissue_df$norm_ref_count, by = list(gene = tissue_df$gene), FUN = sum)
colnames(tissue_df3)[2] <- "ref_allele_count"

# Merge data frames and filter out genes with very low counts
tissue_df_list <- list(tissue_df2, tissue_df3, tissue_df1)
tissue_df_final <- Reduce(function(x, y) merge(x, y, all = TRUE), tissue_df_list, accumulate = FALSE)
tissue_df_final <- tissue_df_final %>% filter(total_count > 3)

# Round counts and calculate allele ratios
tissue_df_final$alt_allele_count <- round(tissue_df_final$alt_allele_count, digits = 0)
tissue_df_final$ref_allele_count <- round(tissue_df_final$ref_allele_count, digits = 0)
tissue_df_final$total_count <- tissue_df_final$alt_allele_count + tissue_df_final$ref_allele_count

tissue_df_final$inv_allele_ratio <- tissue_df_final$alt_allele_count / tissue_df_final$total_count
tissue_df_final$inv_allele_ratio <- round(tissue_df_final$inv_allele_ratio, digits = 3)
tissue_df_final$ref_allele_ratio <- tissue_df_final$ref_allele_count / tissue_df_final$total_count
tissue_df_final$ref_allele_ratio <- round(tissue_df_final$ref_allele_ratio, digits = 3)

# Add tissue column
tissue_df_final <- add_column(tissue_df_final, tissue = "Hypothalamus", .after = "gene")
colnames(tissue_df_final)[2] <- "Tissue"

# Perform binomial test for each gene and calculate p-values
bt <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative = "two.sided", conf.level = 0.95)$p.value}
tissue_df_final$pvalue <- mapply(bt, tissue_df_final$alt_allele_count, tissue_df_final$total_count)

# Adjust p-values for multiple testing using the Benjamini-Hochberg method
tissue_df_final$padjust <- p.adjust(tissue_df_final$pvalue, method = "BH", n = length(tissue_df_final$pvalue))

# Save results to CSV file
write.csv(tissue_df_final, file = "allele_specific_expression/ASE_results/tissues/faeder/results_tissue_faeder_Hyp_ASE.csv", 
          row.names = FALSE)
