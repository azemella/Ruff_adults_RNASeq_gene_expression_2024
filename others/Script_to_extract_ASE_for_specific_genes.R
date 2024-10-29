# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression/allele_specific_expression_analysis/ASE_results/individuals/satellite/")
dir <- "/home/alex-zemella/Documents/ruff_adults_gene_expression/allele_specific_expression_analysis/ASE_results/individuals/satellite/"

# Load required libraries
library("tidyverse")

# Function to find inversion gene of interest
file_list <- list.files(dir, full.names = TRUE, pattern = ".csv")
result_df <- data.frame()

for (file_path in file_list) {
  input_file <- read.csv(file_path, row.names = 1)
  filename <- tools::file_path_sans_ext(basename(file_path))
  input_file$alt_allele_allelic_ratio <- round(input_file$alt_allele_allelic_ratio, digits = 2)
  row_index <- which(rownames(input_file) == "HSD17B2")
  if (length(row_index) > 0) {
    hsd17b2_value <- input_file[row_index, "alt_allele_allelic_ratio"]
    hsd17b2_df <- data.frame(filename = filename, alt_allele_allelic_ratio = hsd17b2_value)
    result_df <- rbind(result_df, hsd17b2_df)
  }
}

result_df$filename <- sub("_ASEReadCounter_results_ASE$", "", result_df$filename)
colnames(result_df)[1] <- "Sample"
result_df$Morph <- "Satellite"
result_df$Tissue <- sub(".*_([^_]+)_\\d+$", "\\1", result_df$Sample)
colnames(result_df)[2] <- "HSD17B2"

result_df <- result_df %>%
  mutate(Tissue = recode(Tissue,
                         "GON" = "Gon",
                         "ADR" = "Adr",
                         "PIT" = "Pit",
                         "TNA" = "A+TnA",
                         "VTA" = "VTA+SN",
                         "RAP" = "Rap",
                         "LIV" = "Liv",
                         "HYP" = "Hyp",
                         "POA" = "POA",
                         "LS" = "LS",
                         "NC" = "NC"))

# Save results as csv
write.csv(result_df, file = "/home/alex-zemella/Documents/ruff_adults_gene_expression/others/HSD17B2_ASE_Satellite.csv", row.names = FALSE)
