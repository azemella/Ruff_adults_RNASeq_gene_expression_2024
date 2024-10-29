# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyverse")
library("PMCMRplus")

# Read hormone data
data <- read.csv(file = "data/hormones_and_behaviors_data/Testes_and_blood_hormone_data.csv")

# Calculate ratios
data <- data %>% mutate(gon_A4overT = (gon_A4/gon_T))
data <- data %>% mutate(blo_A4overT = (blo_A4/blo_T))

# Normality test
shapiro.test(data$gon_A4overT)
# Shapiro-Wilk normality test
# data:  data$gon_A4overT
# W = 0.67459, p-value = 1.155e-07

# Normality test
shapiro.test(data$blo_A4overT)
# Shapiro-Wilk normality test
# data:  data$blo_A4overT
# W = 0.84107, p-value = 0.001198

# Kruskal-wallis test
kruskal.test(data$gon_A4overT ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$gon_A4overT by data$morph
# Kruskal-Wallis chi-squared = 20.968, df = 2, p-value = 2.799e-05

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$gon_A4overT, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$gon_A4overT and data$morph 
# Fae     Ind    
# Ind 0.00026 -      
#   Sat 1.00000 3.2e-05
# P value adjustment method: bonferroni 

# Kruskal-wallis test
kruskal.test(data$blo_A4overT ~ data$morph)
# Kruskal-Wallis rank sum test
# data:  data$blo_A4overT by data$morph
# Kruskal-Wallis chi-squared = 17.308, df = 2, p-value = 0.0001745

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(data$blo_A4overT, data$morph, p.adjust.method = "bonferroni")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  data$blo_A4overT and data$morph 
# Fae     Ind    
# Ind 0.00031 -      
#   Sat 1.00000 0.00014
# P value adjustment method: bonferroni 
