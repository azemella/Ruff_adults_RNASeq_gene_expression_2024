# Authors: Jasmine L. Loveland (jasmine.loveland@univie.ac.at), Alex Zemella (alex.zemella@bi.mpg.de)
# Affiliations: University of Vienna (JL), Max Planck Institute for Biological Intelligence (AZ)

# Set the working directory
setwd("/home/alex-zemella/Documents/ruff_adults_gene_expression")

# Load required libraries
library("tidyverse")
library("PMCMRplus")

# Read data files for aggression and courtship morphs
agg <- read.table("data/hormones_and_behaviors_data/Aggression_morphs_RNAseq.txt", header = TRUE, sep = "\t")
court <- read.table("data/hormones_and_behaviors_data/Courtship_morphs_RNAseq.txt", header = TRUE, sep = "\t")

# Normality test
shapiro.test(agg$aggressionrate)
# Shapiro-Wilk normality test
# data:  agg$aggressionrate
# W = 0.6803, p-value = 2.882e-06

# Normality test
shapiro.test(court$courtshiprate)
# Shapiro-Wilk normality test
# data:  court$courtshiprate
# W = 0.87853, p-value = 0.00536

# Kruskal-wallis test
kruskal.test(agg$aggressionrate ~ agg$morph)
# Kruskal-Wallis rank sum test
# data:  agg$aggressionrate by agg$morph
# Kruskal-Wallis chi-squared = 21.03, df = 2, p-value = 2.712e-05

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(agg$aggressionrate, agg$morph, p.adjust.method = "BH")
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  agg$aggressionrate and agg$morph 
# Fae    Ind   
# Ind 0.0006 -     
#   Sat 0.9315 0.0006
# P value adjustment method: BH 

# Kruskal-wallis test
kruskal.test(court$courtshiprate ~ court$morph)
# Kruskal-Wallis rank sum test
# data:  court$courtshiprate by court$morph
# Kruskal-Wallis chi-squared = 16.983, df = 2, p-value = 0.0002052

# Pairwise comparisons with p-value adjustment
pairwise.wilcox.test(court$courtshiprate, court$morph, p.adjust.method = "BH")
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# data:  court$courtshiprate and court$morph 
# Fae    Ind   
# Ind 0.0028 -     
#   Sat 0.0025 0.0274
# P value adjustment method: BH 
