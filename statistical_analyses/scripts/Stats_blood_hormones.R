# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyverse")
library("PMCMRplus")

# Read hormone data
data <- read.csv(file = "data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv")

# Remove rows with NAs
data <- data[!is.na(data$blo_T), ]

# Normality test
shapiro.test(data$blo_T)
# Shapiro-Wilk normality test
# data:  data$blo_T
# W = 0.4617, p-value = 2.408e-10

# Normality test
shapiro.test(data$blo_A4)
# Shapiro-Wilk normality test
# data:  data$blo_A4
# W = 0.83748, p-value = 0.001028

# Kruskal-wallis test
kruskal.test(data$blo_T ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$blo_T by data$morph
# Kruskal-Wallis chi-squared = 14.384, df = 2, p-value = 0.0007524

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$blo_T, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$blo_T and data$morph 
# Fae    Ind   
# Ind 0.0016 -     
#   Sat 1.0000 0.0065
# P value adjustment method: bonferroni 

# Kruskal-wallis test
kruskal.test(data$blo_A4 ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$blo_A4 by data$morph
# Kruskal-Wallis chi-squared = 17.352, df = 2, p-value = 0.0001706

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$blo_A4, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$blo_A4 and data$morph 
# Fae     Ind    
# Ind 0.00031 -      
#   Sat 1.00000 0.00014
# P value adjustment method: bonferroni 
