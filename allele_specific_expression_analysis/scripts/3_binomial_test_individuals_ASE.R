# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set working directory and define directory path
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required libraries
library("tidyverse")
library("magrittr")
library("extraoperators")

# Run this script for "faeder" and "satellite" files separately
# Replace the inversion morph name in lines #17

# Get a list of files in the specified directory
my_files <- list.files(path = "allele_specific_expression/Edited_GATK_ASEReadCounter_output/faeder/", 
                       pattern = ".csv", full.names = TRUE, recursive = FALSE)

# Define a function to measure allele-specific expression
measure_ase <- function(one_file) {
  # Read data from the CSV file
  data_all <- read.csv(one_file, header = TRUE)

  # Aggregate total count, alternate allele count, and reference allele count by gene
  data_all1 <- aggregate(data_all$total_count, by = list(gene = data_all$gene), FUN = sum)
  colnames(data_all1)[2] <- "total_count"

  data_all2 <- aggregate(data_all$alt_count, by = list(gene = data_all$gene), FUN = sum)
  colnames(data_all2)[2] <- "alt_allele_count"

  data_all3 <- aggregate(data_all$ref_count, by = list(gene = data_all$gene), FUN = sum)
  colnames(data_all3)[2] <- "ref_allele_count"

  # Combine aggregated data frames
  df_list <- list(data_all2, data_all3, data_all1)
  data_all_final <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
  
  # Calculate allelic ratio
  data_all_final$alt_allele_allelic_ratio <- data_all_final$alt_allele_count / data_all_final$total_count
  
  # Perform binomial test for each gene and calculate p-values
  bt <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative = "two.sided", conf.level = 0.95)$p.value}
  data_all_final$pvalue <- mapply(bt, data_all_final$alt_allele_count, data_all_final$total_count)
  
  # Adjust p-values for multiple testing using the Benjamini-Hochberg method
  data_all_final$padjust <- p.adjust(data_all_final$pvalue, method = "BH", n = length(data_all_final$pvalue))
  
  # Write the results to a new CSV file
  write.csv(data_all_final, 
            file = sub("\\_edited.csv", "\\_results_ASE.csv", one_file),
            row.names = FALSE)
  
  # Return the final data frame
  return(data_all_final)
}

# Apply the measure_ase function to each file in the list of files
res <- lapply(my_files, function(i) measure_ase(i)) 
