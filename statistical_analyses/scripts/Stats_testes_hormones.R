# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyverse")
library("PMCMRplus")

# Read hormone data
data <- read.csv(file = "data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv")

# Normality test
shapiro.test(data$gon_T)
# Shapiro-Wilk normality test
# data:  data$gon_T
# W = 0.85477, p-value = 0.0002455

# Normality test
shapiro.test(data$gon_A4)
# Shapiro-Wilk normality test
# data:  data$gon_A4
# W = 0.83895, p-value = 0.000107

# Kruskal-wallis test
kruskal.test(data$gon_T ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$gon_T by data$morph
# Kruskal-Wallis chi-squared = 12.311, df = 2, p-value = 0.002122

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$gon_T, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$gon_T and data$morph 
# Fae    Ind   
# Ind 0.0069 -     
#   Sat 1.0000 0.0116
# P value adjustment method: bonferroni 

# Kruskal-wallis test
kruskal.test(data$gon_A4 ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$gon_A4 by data$morph
# Kruskal-Wallis chi-squared = 26.075, df = 2, p-value = 2.178e-06

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$gon_A4, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$gon_A4 and data$morph 
# Fae     Ind    
# Ind 1.9e-06 -      
#   Sat 1       1.4e-06
# P value adjustment method: bonferroni 