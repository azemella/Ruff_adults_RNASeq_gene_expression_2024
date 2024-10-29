# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression"

# Load required library
library("tidyverse")

### Run this script for "faeder" and "satellite" files separately
### Replace the inversion morph name in lines #57 and #62

# Define a function to process each file
process_file <- function(file_path) {
  
  # Read the input file
  input_file <- read.csv(file_path)
  file_name <- basename(file_path)
  
  # Read library size data
  lib_size <- read.csv("metadata/metadata.csv")
  lib_size <- lib_size %>% distinct(Sample_ID, .keep_all = TRUE)
  lib_size <- lib_size[, c("Sample_ID", "Lib_size")]
  lib_size <- lib_size %>% remove_rownames %>% column_to_rownames(var = "Sample_ID")
  
  # Update row names to match the file names
  rownames(lib_size) <- paste(rownames(lib_size), "_ASEReadCounter_edited.csv", sep = "")
  
  # Find matching library size for the input file
  matching_value <- rownames(lib_size)[grepl(file_name, rownames(lib_size))]
  corresponding_value <- lib_size[matching_value, ]
  
  # Add library size column to input file
  input_file$library_size <- corresponding_value
  
  # Reorder columns
  input_file <- relocate(input_file, library_size, .before = ref_count)
  
  # Normalize ref_count and alt_count
  input_file$norm_ref_count <- input_file$ref_count / (input_file$library_size / 1e7)
  input_file$norm_ref_count <- round(input_file$norm_ref_count, digits = 2)
  input_file$norm_alt_count <- input_file$alt_count / (input_file$library_size / 1e7)
  input_file$norm_alt_count <- round(input_file$norm_alt_count, digits = 2)
  
  # Calculate total_count
  input_file$total_count <- input_file$norm_ref_count + input_file$norm_alt_count
  input_file$total_count <- round(input_file$total_count, digits = 2)
  
  # Remove unnecessary columns
  input_file <- input_file[,-c(5,7)]
  
  # Reorder columns
  input_file <- input_file %>% relocate(norm_ref_count, .before = total_count ) %>% relocate(norm_alt_count, .before = total_count )
  
  # Save the modified data to a new file (optional)
  output_file_path <- file.path("allele_specific_expression_analysis/Normalization/faeder/", file_name)
  write.csv(input_file, file = output_file_path, row.names = FALSE)
}

# Run this code on the shell command afterwards
# find /home/alex-zemella/Documents/PhD/results/allele_specific_expression/Edited_GATK_ASEReadCounter_output/faeder/ -type f -name '*.csv' -exec Rscript -e "source('2-normalize_ASE_counts.R'); process_file('{}')" \;
